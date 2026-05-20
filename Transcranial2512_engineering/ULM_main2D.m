clear; clc; close all;

%% 选择数据并设置保存路径
projectName = 'Transcranial2512';
expName     = '20250411flow';
baseDir     = fileparts(mfilename('fullpath'));

% 路径和输出目录创建统一交给 setupULMPaths
[dataDir,localPath,trackspath,savingpath,fileList,allFileNames,nBuffers] = setupULMPaths(baseDir, expName);

%% 参数设置
tempFile = fullfile(dataDir, allFileNames{1}); 
loadedData = load(tempFile);
fieldNames = fieldnames(loadedData);
IQ = loadedData.IQData; % 需要先确认变量名
% tempData = struct2cell(load(targetFile)); 
% IQ = tempData{1}; 

[nAxial, nLateral, nFrames] = size(IQ);clear IQ;
pitch        = 0.2/1000; % 阵元间距
frameRate    = 1000; % 帧率 (Hz)
nElements    = 128; % 阵元数量
centerFreq   = 15e6; % 中心频率 (15 MHz)
cValue       = 1540; % 声速
lambda       = cValue / centerFreq; % 波长
resFactor    = 10; % 超分辨放大倍数

ULM = struct();
ULM.size                   = [nAxial, nLateral, nFrames];
ULM.numberOfParticles      = 200;
ULM.res                    = resFactor;% 超分辨因子。最终生成的超分辨图像像素比原始 IQ 图像精细 10 倍。
ULM.interpFactor           = 1 / resFactor;% 定位精度因子。在子像素定位时，将原始像素切分的粒度。
ULM.svdCutoff              = round([10, nFrames]);% round四舍五入到最接近的整数
ULM.max_linking_distance   = 2;% 最大连接距离
ULM.min_length             = 7;% 最小轨迹长度
ULM.fwhm                   = [3, 3];% 半高全宽 (PSF)
ULM.max_gap_closing        = 0;% 跳帧补偿,表示只要断开一帧，轨迹立即终止。
ULM.scale                  = [1, 1, 1/frameRate]; % [z, x, t] 的转换系数
ULM.numberOfFrameProcessed = nFrames;% 需要实际处理多少帧

% 滤波器配置 (Butterworth)
% ULM.butter.samplingFreq / 2 计算的是奈奎斯特频率（Nyquist Frequency），最高只能准确捕捉到的信号变化。
ULM.butter.cutoffFreq   = [50, 249];
ULM.butter.samplingFreq = frameRate;
[butB, butA] = butter(2, ULM.butter.cutoffFreq / (ULM.butter.samplingFreq/2), 'bandpass'); % 即只让 50-249Hz 之间的信号通过
ULM.parameters.NLocalMax = 4; % 在识别微泡中心时，算法会以一个像素为中心，检查周围 4x4 像素范围内，它是否是能量最强的点

PData.Size   = [nAxial, nLateral, nFrames];
PData.PDelta = [1, 1, 1]; % [dx, dy, dz]，步长
PData.Origin = [0, PData.Size(2)/2*PData.PDelta(2), 0];% 原点位置，将侧向轴 (x) 的中心设为 0，符合探头物理排列

% 计算物理尺寸 (单位: mm)
widthM  = floor(nLateral) * lambda; % floor函数向下取整
depthM  = floor(nAxial) * lambda;
xAxis_mm = linspace(0, widthM, nLateral) .* 1000;
zAxis_mm = linspace(0, depthM, nAxial) .* 1000;
%% Location and Tracking
% for hhh = 1 : nBuffers 
for hhh = 1 : 1
% parfor hhh = 1:min(999,nBuffers)
    fprintf('Processing block %d/%d: %s\n', hhh, nBuffers, allFileNames{hhh});
    dataLoad = load(fullfile(dataDir, allFileNames{hhh}));
    temp = dataLoad.IQData; % 需要先确认变量名
    temp([1, end], :, :) = 0;
    temp(:, [1, end], :) = 0;% 轴向边缘清零,第一/最后一行列清零，图像边缘常会产生伪影
    
    IQ_filt = SVDfilter(temp, ULM.svdCutoff);
%     IQ_filt = filter(butB,butA,IQ_filt,[],3);
    IQ_filt(~isfinite(IQ_filt))=0; % 检查矩阵中的每个数是否都是有限值
    
    motionCorrection.delta = zeros(800, 2); % 运动补偿为0
    
    trackFile = fullfile(trackspath, sprintf('Tracks_%03d.mat', hhh));% sprintf格式化数据并将其返回为字符串
    localFile = fullfile(localPath, sprintf('Locals_%03d.mat', hhh)); 
    multiULM(IQ_filt,ULM,PData,motionCorrection.delta, ...
                        'savingTrackfilename', trackFile, ...
                        'savingLocalfilename', localFile);
    
    clear temp IQ_filt dataLoad;
end

MatOut= [];% 密度图
MatOut_vel = MatOut;% 速度图           
MatOut_z = [];% z向速度
MatOut_x = [];% x向速度
index = 1;
MatOutBubble = [];
MatOutSat = [];

% parfor hhh = 1:min(999,Nbuffers)
% for hhh = 1:nBuffers %1:Nbuffers 
for hhh = 1:1
    trackLength = [];   
    trackFullFile = fullfile(trackspath, sprintf('Tracks_%03d.mat', hhh));
    localFullFile = fullfile(localPath, sprintf('Locals_%03d.mat', hhh));
    
    % 检查文件是否存在，避免报错中断
    if ~exist(trackFullFile, 'file')
        fprintf('警告: 找不到文件 %s，跳过...\n', trackFullFile);
        continue;
    end
    % 加载文件
    load(trackFullFile, 'Track_raw','Track_interp','ProTime');
    load(localFullFile);
    
    aa = -PData.Origin([3 1])+[1 1]*1; %[1 1 0]
    bb = 1./PData.PDelta([3 1]);%[1 1 1 1 1 1]
    aa(3) = 0;
    bb(3:6) = 1;
    
%     检查这个大矩阵的列数
    if size(cell2mat(Track_raw),2) == 4
        MatOutBubble(index,1) = 0;
        index = index + 1;
        continue;
    end
    
%     [Track_raw,Track_interp] = kalman_filter(Track_raw,ULM); 
%     [Track_raw,Track_interp] = angle_constrain(Track_raw,Track_interp);
%     [Track_raw,Track_interp] = accela_constrain(Track_raw,Track_interp);
    
    Track_matout = Track_interp;
    Track_matout = cellfun(@(x) (x(:,[1 2 3 4 5 7]).* bb),Track_matout,'UniformOutput',0);
    Track_raw = cellfun(@(x) (x(:,[1 2 3 4 5 7]).* bb),Track_raw,'UniformOutput',0);
    for i = 1:length(Track_raw)
        singleTrack = Track_raw{i};
        trackLength(i) = length(singleTrack);% 统计每一条轨迹到底活了多少帧
    end
    
    meanLength(hhh) = mean(trackLength);
    trackNumber(hhh) = length(Track_raw);
    trackBubble(hhh) = sum(trackLength(:))/size(MatTracking,1);
    fprintf('Block %03d | 轨迹数: %4d | 平均长度: %5.2f 帧 | 利用率: %5.2f%%\n', ...
        hhh, trackNumber(hhh), meanLength(hhh), trackBubble(hhh)*100);
    
    [MatOut_i,MatOut_vel_i,MatOut_z_i,MatOut_x_i, MatOut_int_i] = ULM_Track2MatOut(Track_matout,ULM.res*[PData(1).Size(1) PData(1).Size(2)]+[1 1]*1,'mode','2D_vel_z','2_velmean',ULM);
    s_temp_size = size(MatOut);
    if s_temp_size == [0,0]
        MatOut = zeros(size(MatOut_i));
        MatOut_vel = zeros(size(MatOut_vel_i));
        MatOut_z = zeros(size(MatOut_z_i));
        MatOut_x = zeros(size(MatOut_x_i));
        MatOut_int = zeros(size(MatOut_int_i));
        clear s_temp_size
    end
    clear Track_matout
    
    MatOut_vel = MatOut_vel.*MatOut+MatOut_vel_i.*MatOut_i; % weighted summation
    MatOut_z = MatOut_z.*MatOut+MatOut_z_i.*MatOut_i;
    MatOut_x = MatOut_x.*MatOut+MatOut_x_i.*MatOut_i;
    MatOut_int = MatOut_int.*MatOut + MatOut_int_i.*MatOut_i;
    
    MatOut = MatOut+MatOut_i;
    MatOut_vel(MatOut>0) = MatOut_vel(MatOut>0)./MatOut(MatOut>0); % average velocity
    MatOut_z(MatOut>0) = MatOut_z(MatOut>0)./MatOut(MatOut>0);
    MatOut_x(MatOut>0) = MatOut_x(MatOut>0)./MatOut(MatOut>0);
    MatOut_int(MatOut>0) = MatOut_int(MatOut>0)./MatOut(MatOut>0);
    
    MatOutBubble(index,1) = sum(MatOut_i(:));
%     MatOutBubble(index,1) = nnz(MatOut_i(:));
    MatOutSat(index,1) = nnz(MatOut>0); 
%     MatOut
    index = index + 1;
end
    
    UF.FrameRateUF = 1000;
    UF.F0 = 15e6;
clear  Track_interp Track_count Track_matout 
save(fullfile(savingpath, 'MatOut_multi.mat'), 'MatOut', 'MatOut_vel', 'MatOutSat', 'ULM', 'PData', 'UF');

% trackLengthKSVD(k) = mean(meanLength);
% bubbleTrack(k) = sum(MatOut(:));
% trackedPercent(k) = mean(trackBubble);
% MatOutSatK{k} = MatOutSat;

% end

% load([savingpath 'MatOut_multi_slice0mm']);
c = 1540;
f0 = 15e6;
lambda = c/f0;
Nelements = 128;

[M,N] = size(MatOut);
Width = lambda * floor(N/10);
Depth = lambda * floor(M/10);
z = linspace(0,Depth,M).*1000;% 探头尺寸
x = linspace(0,Width,N).*1000;
% MatZoom = [373 590 735 930];

figure(2)
IntPower = 1/2;
SigmaGauss=0.3;
im=imagesc(x,z,imgaussfilt(MatOut.^IntPower,.01));axis image
if SigmaGauss>0,im.CData = imgaussfilt(im.CData,SigmaGauss);end % 对图像进行高斯模糊

title('ULM intensity display')
colormap(gca,hot(128))
% clbar = colorbar;
caxis(caxis*.8)  % 饱和度截断
% clbar.Label.String = 'number of counts';
% clbar.TickLabels = round(clbar.Ticks.^(1/IntPower),1);
% xlabel('\lambda');ylabel('\lambda') 
ca = gca;ca.Position = [.05 .05 .8 .9];

% im.CData = min(im.CData,10);caxis([0 10]) 
axis image
ylabel('Axial Distance(mm)')
xlabel('Lateral Distance(mm)')
% axis offSR = imgaussfilt(MatOut.^IntPower,.01);


% MatOut_vel = MatOut_vel./10;
% figure(3)
% vmax_disp  = ceil(quantile(MatOut_vel(abs(MatOut_vel)>0),.98)/10)*10;% 找到最红颜色代表的速度
% IntPower = 1/2;
% lambda = 1540/15e6 * 1e3;
% ULM.SRscale = 10;
% 
% clf,set(gcf,'Position',[652 393 941 585]);
% clbsize = [180,50];
% Mvel_rgb = MatOut_vel/vmax_disp; % normalization
% Mvel_rgb(1:clbsize(1),1:clbsize(2)) = repmat(linspace(1,0,clbsize(1))',1,clbsize(2)); % 180 x 50 的矩形块，直接覆盖到图像数据
% Mvel_rgb = Mvel_rgb.^(1/2.5);%伽马校正，非线性地提亮图像
% Mvel_rgb(Mvel_rgb>1)=1;% 所有大于 1 的值设为 1
% Mvel_rgb = abs(Mvel_rgb);% 确保没有负数
% Mvel_rgb = imgaussfilt(Mvel_rgb,.7);% 使用高斯滤波对矩阵进行轻微模糊处理
% Mvel_rgb = ind2rgb(round(Mvel_rgb*256),jet(256)); % convert ind into RGB
% 
% MatShadow = MatOut;
% MatShadow = MatShadow./max(MatShadow(:)*.3);% 把图像亮度放大约 3.3 倍
% MatShadow(MatShadow>1)=1;
% MatShadow(1:clbsize(1),1:clbsize(2))=repmat(linspace(0,1,clbsize(2)),clbsize(1),1);
% Mvel_rgb = Mvel_rgb.*(MatShadow.^IntPower);% 融合
% Mvel_rgb = brighten(Mvel_rgb,.2);% 非线性地提亮整个图像的颜色映射
% BarWidth = round(1./(ULM.SRscale*lambda)); % 1 mm在图像中对应多少个像素
% Mvel_rgb(size(MatOut,1)-50+[0:3],60+[0:BarWidth],1:3)=1;% 手动绘制比例尺
% imshow(Mvel_rgb);axis on
% title(['Velocity magnitude (0-' num2str(vmax_disp) 'mm/s)'])
% ca = gca;ca.Position = [.05 .05 .8 .9];
% 
% figure(4)
% MatOut_zdir = MatOut_vel;
% velColormap = cat(1,flip(flip(hot(256),1),2),hot(256)); % custom velocity colormap
% velColormap = velColormap(5:end-5,:); % remove white parts红正蓝负
% IntPower = 1/2;
% im=imagesc(x,z,(MatOut).^IntPower.*sign(imgaussfilt(MatOut_zdir,.3)));%(坐标，亮度，符号方向）
% im.CData = im.CData - sign(im.CData)/2;axis image
% title(['ULM intensity display with axial flow direction'])
% colormap(gca,velColormap)
% caxis([-1 1]*max(caxis)*.8) % add saturation in image
% clbar = colorbar;clbar.Label.String = 'Count intensity';
% ca = gca;ca.Position = [.05 .05 .8 .9];
% figure(6)
% vmax_disp  = ceil(quantile(MatOut_int(abs(MatOut_int)>0),.98)/10)*10;
% IntPower = 1/2;
% lambda = 1;
% ULM.SRscale = 10;
% clf,set(gcf,'Position',[652 393 941 585]);
% clbsize = [180,50];
% Mvel_rgb = MatOut_int/vmax_disp; % normalization
% Mvel_rgb(1:clbsize(1),1:clbsize(2)) = repmat(linspace(1,0,clbsize(1))',1,clbsize(2)); % add velocity colorbar
% Mvel_rgb = Mvel_rgb.^(1/1.5);Mvel_rgb(Mvel_rgb>1)=1;
% Mvel_rgb = abs(Mvel_rgb);
% Mvel_rgb = imgaussfilt(Mvel_rgb,.5);
% Mvel_rgb = ind2rgb(round(Mvel_rgb*256),hot(128)); % convert ind into RGB
% 
% MatShadow = MatOut;
% MatShadow = MatShadow./max(MatShadow(:)*.4);
% MatShadow(MatShadow>1)=1;
% MatShadow(1:clbsize(1),1:clbsize(2))=repmat(linspace(0,1,clbsize(2)),clbsize(1),1);
% Mvel_rgb = Mvel_rgb.*(MatShadow.^IntPower);
% Mvel_rgb = brighten(Mvel_rgb,.3);
% % BarWidth = round(1./(ULM.SRscale*lambda)); % 1 mm
% % Mvel_rgb(size(MatOut,1)-50+[0:3],60+[0:BarWidth],1:3)=1;
% imshow(Mvel_rgb);axis on
%  
% 


% BarWidth = round(1./(ULM.scale(2)*lambda)); % 1 mm
% im.CData(size(im.CData,1)-2,3+[0:BarWidth])=max(caxis);


% figure(7)
% % t = linspace(1,round(length(MatOutSat)/20),length(MatOutBubble));
% % plot(t, MatOutSat)
% % plot(smooth(MatOutSat,15));
% % 修改后（使用官方内置的平滑函数）
% plot(smoothdata(MatOutSat, 'movmean', 15));
% xlabel('Data Batch');
% ylabel('Number of Bubbles Tracked')


% 
% MatOut_structure = MatOut;
% MatOut_velocity = cat(3,MatOut_z,MatOut_x);

