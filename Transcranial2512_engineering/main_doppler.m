clear; clc; close all;

projectName = 'Transcranial2512';
expName     = '20250411flow';
baseDir     = fileparts(mfilename('fullpath'));
% 路径和输出目录创建统一交给 setupULMPaths，避免诊断脚本重复维护路径逻辑。
[dataDir,localPath,trackspath,savingpath,fileList,allFileNames,nBuffers] = setupULMPaths(baseDir, expName);
% T = struct2table(fileList);
% disp(T);

targetFile = fullfile(dataDir, allFileNames{1}); 
loadedData = load(targetFile);
fieldNames = fieldnames(loadedData);
% tempData = struct2cell(load(targetFile)); 
% IQ = tempData{1}; 
% clear tempData;% 清理中间大变量，释放内存
if isfield(loadedData, 'IQData')
    IQ = loadedData.IQData;
    fprintf('Detected input variable: IQData\n');
elseif isfield(loadedData, 'IQ')
    IQ = loadedData.IQ;
    fprintf('Detected input variable: IQ\n');
elseif ~isempty(fieldNames)
    % 备选方案：抓取文件中的第一个变量
    IQ = loadedData.(fieldNames{1});
    warning('Standard input variable name not found. Using first variable: %s', fieldNames{1});
else
    error('Input file is empty or does not contain a valid variable: %s', allFileNames{1});
end
[nAxial, nLateral, nFrames] = size(IQ);
pitch        = 0.2/1000; % 阵元间距
frameRate    = 1000; % 帧率 (Hz)
centerFreq   = 15e6; % 中心频率 (15 MHz)
cValue       = 1540; % 声速
lambda       = cValue / centerFreq; % 波长
nElements    = 128; % 阵元数量

maxBubbles      = 200;
resFactor       = 10;
connectDist     = 2;
minTrackLen     = 7;
bandpassFreq    = [50, 249]; 
nLocalMax       = 4;

ULM = struct();
ULM.size                   = [nAxial, nLateral, nFrames];
ULM.numberOfParticles      = maxBubbles;
ULM.res                    = resFactor;% 超分辨因子。最终生成的超分辨图像像素比原始 IQ 图像精细 10 倍。
ULM.interpFactor           = 1 / resFactor;% 定位精度因子。在子像素定位时，将原始像素切分的粒度。
ULM.svdCutoff              = round([10, nFrames]);% round四舍五入到最接近的整数
ULM.max_linking_distance   = connectDist;% 最大连接距离
ULM.min_length             = minTrackLen;% 最小轨迹长度
ULM.fwhm                   = [3, 3];% 半高全宽 (PSF)
ULM.max_gap_closing        = 0;% 跳帧补偿,表示只要断开一帧，轨迹立即终止。
ULM.scale                  = [1, 1, 1/frameRate]; % [z, x, t] 的转换系数
ULM.numberOfFrameProcessed = nFrames;% 需要实际处理多少帧

% 滤波器配置 (Butterworth)
% ULM.butter.samplingFreq / 2 计算的是奈奎斯特频率（Nyquist Frequency），最高只能准确捕捉到的信号变化。
ULM.butter.cutoffFreq   = bandpassFreq;
ULM.butter.samplingFreq = frameRate;
[butB, butA] = butter(2, ULM.butter.cutoffFreq / (ULM.butter.samplingFreq/2), 'bandpass'); % 即只让 50-249Hz 之间的信号通过
ULM.parameters.NLocalMax = nLocalMax; % 在识别微泡中心时，算法会以一个像素为中心，检查周围 4x4 像素范围内，它是否是能量最强的点

PData.Size   = [nAxial, nLateral, nFrames];
PData.PDelta = [1, 0, 1]; % [x, y, z]，步长
PData.Origin = [0, PData.Size(2)/2*PData.PDelta(2), 0];% 将侧向轴 (x) 的中心设为 0，符合探头物理排列

% 计算物理尺寸 (单位: m)
widthM  = floor(nLateral) * lambda; % floor函数向下取整
depthM  = floor(nAxial) * lambda;
xAxis_mm = linspace(0, widthM, nLateral) .* 1000;
zAxis_mm = linspace(0, depthM, nAxial) .* 1000;
IQ_filt = SVDfilter(IQ,[10, 400]);
% IQ_filt = filter(butB,butA,IQ_filt,[],3);
IQ_filt(~isfinite(IQ_filt))=0;% 检查矩阵里的每一个数字是不是一个"有限的正常数字",NaN/Inf

BullesSNR = abs(IQ_filt(:,:,200));
LocalMax = imregionalmax(BullesSNR);
ValMax = sort(BullesSNR(LocalMax),'descend');% BullesSNR(LocalMax)自动拉平成一维向量
noise = mean(BullesSNR(:));% 使用冒号操作符将矩阵展开成一个一维列向量
thresh = noise * (10^(9/20));% 规定只有信号强度比平均背景噪声高出 9 dB（分贝）的亮点，才算作真正的微泡
mask = (ValMax >= thresh);

% ULM.numberOfParticles = round(mean(numberOfParticle(:)));
SNRmean = 20 * log10(mean(ValMax(1:10))/mean(BullesSNR(:)));% 平均信噪比
% clear BullesSNR LocalMax ValMax;

figure(1)
dB = 20 * log10(abs(IQ_filt));
dB = dB - max(dB(:));
imagesc(xAxis_mm,zAxis_mm,dB(:,:,200),[-50 0]);
colormap gray;colorbar;
axis equal;axis tight;
xlabel('Axial Distance(mm)')
ylabel('Lateral Distance(mm)')
set(gca,'FontSize',15,'FontWeight','bold','FontName','Times New Roman');% 修改坐标轴
% ========== 优化版：奇异值分布 + 斜率法自动找拐点 ==========
[nx, nz, nt] = size(IQ);
Casorati = reshape(IQ, nx * nz, nt);

% 计算奇异值
s = svd(Casorati, 'econ'); 

% 绘图
figure('Name', 'Singular Value Distribution', 'Color', 'w');
semilogy(1:length(s), s, 'b-', 'LineWidth', 2);
hold on; grid on;

% ===================== 斜率法：自动找拐点 =====================
% 原理：计算相邻奇异值的斜率，斜率突然变小的位置 = 噪声拐点
dy = diff(log(s));  % 对数域斜率（最适合超声数据）
[~, cutoff_idx] = max(abs(diff(abs(dy)))); % 自动找拐点

% 标记自动计算的最佳截断值
xline(cutoff_idx, 'r--', ['Best Cutoff = ', num2str(cutoff_idx)], ...
    'LineWidth', 2, 'LabelOrientation', 'horizontal', ...
    'FontSize', 12, 'FontWeight', 'bold');

% 图像设置（删除冗余代码）
title('Singular Value Energy Distribution', 'FontSize', 14);
xlabel('Singular Value Index', 'FontSize', 12);
ylabel('Magnitude (Log Scale)', 'FontSize', 12);
xlim([1, 300]); % 只看前60个足够，后面全是噪声
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

% 输出命令行结果
fprintf('Auto-selected upper SVD cutoff = %d\n', cutoff_idx);

% --- 设定起始补偿位置 ---
startDepth_mm = 5; % 假设颅骨大约在 5mm 厚度，从 5mm 后开始补偿
[~, zStart_idx] = min(abs(zAxis_mm - startDepth_mm)); % 找到对应的像素索引

% 重新定义计算区域
zDim_ROI = zlim - zStart_idx + 1;
sectionSize = floor(zDim_ROI / nSections);

% 初始化
gainFactors = zeros(nSections, 1);
zCenters = zeros(nSections, 1);

% --- 重新计算全局最大值（避开颅骨强回声，只在大脑区域找基准） ---
brainData = BullesSNR(zStart_idx:end, :);
globalMax = max(brainData(:)); 

for i = 1:nSections
    % 计算相对于 zStart_idx 的范围
    zRange = (zStart_idx - 1) + ((i-1)*sectionSize + 1 : min(i*sectionSize, zDim_ROI));
    sectionData = BullesSNR(zRange, :);
    
    zCenters(i) = mean(zRange); 
    currentMax = max(sectionData(:));
    
    if currentMax > 0
        gainFactors(i) = globalMax / currentMax;
    else
        gainFactors(i) = NaN; % 标记没泡的区域
    end
end

% 处理无效值
gainFactors = fillmissing(gainFactors, 'linear', 'EndValues', 'nearest');

% --- 构建完整的增益曲线 ---
% 1. 在颅骨区域（1 到 zStart_idx），增益设为 1
% 2. 在大脑区域，进行 pchip 插值
zAxis_idx = (1:zDim)';
gainCurve = ones(zDim, 1); % 默认全 1

% 仅在大脑区域应用插值计算出的增益
brain_idx = zStart_idx:zDim;
gainCurve(brain_idx) = pchip(zCenters, gainFactors, brain_idx);

% gainCurve(gainCurve > 15) = 15; %最大范围

figure(2); % 创建新窗口，不覆盖之前的图像
plot(zAxis_mm, gainCurve, 'Color', [0.85 0.33 0.1], 'LineWidth', 2); % 使用深橘色，线宽设为2
grid on;
xlabel('Depth [mm]');
ylabel('Gain');
title('Intensity Gain Compensation Curve');

xlim([zAxis_mm(1), zAxis_mm(end)]);
ylim([0, max(gainCurve) + 2]); % 纵轴留一点余量

set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
BullesSNR_compensated = BullesSNR .* repmat(gainCurve, 1, size(BullesSNR, 2)); % repmat 将曲线扩展到与图像宽度一致

LocalMax = imregionalmax(BullesSNR_compensated);% 之后再使用补偿后的图像进行检测
ValMax = sort(BullesSNR_compensated(LocalMax),'descend');
noise = mean(BullesSNR_compensated(:));% 使用冒号操作符将矩阵展开成一个一维列向量
thresh = noise * (10^(9/20));
mask = (ValMax >= thresh);

% ULM.numberOfParticles = round(mean(numberOfParticle(:)));
SNRmean = 20 * log10(mean(ValMax(1:10))/mean(BullesSNR_compensated(:)));
% clear BullesSNR LocalMax ValMax;

figure(3)
dB = 20 * log10(abs(BullesSNR_compensated));
dB = dB - max(dB(:));
imagesc(xAxis_mm,zAxis_mm,dB,[-40 0]),colormap gray;
colorbar;axis equal;axis tight;
xlabel('Axias Distance(mm)')
ylabel('Lateral Distance(mm)')
set(gca,'FontSize',15,'FontWeight','bold');
set(gca,'FontName','Times New Roman');


figure(5)
PowDop = [];
for hhh=2:2
    path = fullfile(dataDir, allFileNames{hhh});
    tmp = load(path);
    IQ_filt = SVDfilter(tmp.IQData,ULM.svdCutoff);tmp = [];
    IQ_filt = filter(butB,butA,IQ_filt,[],3);
    IQ_filt(~isfinite(IQ_filt))=0;
    PowDop(:,:,end) = sqrt(sum(abs(IQ_filt).^2,3));

    pause(0.1)
end
im=imagesc(mean(PowDop,3).^(1/2));
axis image, colormap(gca,hot(128)),title(['Power Doppler'])
clbar = colorbar;
caxis([10 max(im.CData(:))*.9]);

xlabel('Axial Distance(mm)')
ylabel('Lateral Distance(mm)')
set(gca,'FontSize',15,'FontWeight','bold');
set(gca,'FontName','Times New Roman');
PD = mean(PowDop,3).^(1/2);

