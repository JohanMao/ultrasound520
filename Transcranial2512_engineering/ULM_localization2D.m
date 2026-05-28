function MatTracking = ULM_localization2D(MatIn,ULM,delta)
%ULM_LOCALIZATION2D 此处显示有关此函数的摘要
%   此处显示详细说明
fwhmz = ULM.fwhm(2);
fwhmx = ULM.fwhm(1);

vectfwhmz = -1*round(fwhmz/2):round(fwhmz/2);% [-2, -1, 0, 1, 2]，四舍五入
vectfwhmx = -1*round(fwhmx/2):round(fwhmx/2);
[meshX,meshZ] = meshgrid(vectfwhmx,vectfwhmz); % 生成对应的二维“网格矩阵”
sigGauss_z = vectfwhmz(end)*0+1;%1 微气泡在深度方向上的高斯分布标准差。越小,高斯曲线越尖锐，代表微气泡在图像上是一个非常小的点；越大：高斯曲线越平缓，代表微气泡在图像上是一个模糊的亮团。
sigGauss_x = vectfwhmx(end)*0+1;
Gauss_grid = exp(-(meshZ-0).^2./(2*sigGauss_z^2) - (meshX-0).^2./(2*sigGauss_x^2));% 构建二维高斯分布，中心值最接近 1，越往边缘数值越小

[height,width,numberOfFrames]=size(MatIn);%MatIn = IQ
MatIn = abs(MatIn);% Make sure you work with the intensity matrix
Mat_norm = (MatIn-min(MatIn(:)))./(max(MatIn(:))-min(MatIn(:)));% 离差标准化，确保全图最亮的像素值变为 1
%如果存在异常亮像素，可以先剔除前 0.1% 的极值再进行归一化
info = whos('MatIn');typename = info.class; % 

if ~isfield(ULM,'LocMethod')% 结构体是否有某个“键”
    ULM.LocMethod = 'curvefitting'; % 曲线拟合法
end

if ~isfield(ULM,'parameters') % Create an empty structure for parameters hosting
    ULM.parameters = struct();
end

if strcmp(ULM.LocMethod,'interp')% String Compare
    if ~isfield(ULM.parameters,'InterpMethod')%插值
        ULM.parameters.InterpMethod = 'spline';% 样条插值
    end
    if sum(strcmp(ULM.parameters.InterpMethod,{'bilinear','bicubic'}))
        warning('Faster but pixelated, Weighted Average will be faster and smoother.')% 双线性或双三次插值虽然运算速度快，但结果会呈现“像素化”
    end
end

if ~isfield(ULM.parameters,'NLocalMax')
    if fwhmz==3,ULM.parameters.NLocalMax = 2; %每个小区域只允许有 2 个局部极大值
    else
        ULM.parameters.NLocalMax = 3;
    end
end

%% 根据图像强度选出每一帧中最亮的N个微气泡。
MatInReduced = zeros(height,width,numberOfFrames,typename);
MatInReduced(1+round(fwhmz/2)+1:height-round(fwhmz/2)-1,1+round(fwhmx/2)+1:width-round(fwhmx/2)-1,:) = MatIn(1+round(fwhmz/2)+1:height-round(fwhmz/2)-1, 1+round(fwhmx/2)+1:width-round(fwhmx/2)-1,:);
[height,width,numberOfFrames] = size(MatInReduced);% 复制中心区域，边缘清零

Mat2D = permute(MatInReduced, [1,3,2]); %so that all the frames are in columns 2维和3维互换
Mat2D = reshape(Mat2D,height*numberOfFrames,width); % Concatenate Matrix 3维变2维
mask2D = imregionalmax(Mat2D); clear Mat2D  % Perform imregionalmax 提高效率
mask = reshape(mask2D,height,numberOfFrames,width);clear mask2D % reshape concatenated mask
mask = permute(mask,[1,3,2]); % so that we restore (z,x,t) table
% 候选点矩阵只在局部极大值位置保留强度，其他位置为 0。
IntensityMatrix = MatInReduced.*mask; % 保留局部极大值候选点的强度

% 增益只影响 Top-N 候选点筛选，不直接改变后续亚像素定位使用的原始 ROI 图像。
% 基于 Power Doppler 的皮层 ROI 强度增益补偿。
if isfield(ULM, 'parameters') && isfield(ULM.parameters, 'LocalizationZRange')
    targetZRange = round(ULM.parameters.LocalizationZRange);
else
    targetZRange = [130, 160];
end

% 将目标深度范围限制在当前图像尺寸内，避免越界。
targetZRange(1) = max(1, min(height, targetZRange(1)));
targetZRange(2) = max(targetZRange(1), min(height, targetZRange(2)));

% CandidateZRange 只控制参与 Top-N 竞争的内部区域，避免目标区上下边界抢占定位名额。
candidateZRange = targetZRange;
if isfield(ULM, 'parameters') && isfield(ULM.parameters, 'CandidateZRange')
    candidateZRange = round(ULM.parameters.CandidateZRange);
    candidateZRange = sort(candidateZRange(:).');
    if numel(candidateZRange) ~= 2
        candidateZRange = targetZRange;
    end
end
candidateZRange(1) = max(targetZRange(1), min(targetZRange(2), candidateZRange(1)));
candidateZRange(2) = max(candidateZRange(1), min(targetZRange(2), candidateZRange(2)));

% 在目标区域上下各保留若干行作为缓冲，用于估计增益和避免边界伪影。
processMargin = 5;
processZRange = [ ...
    max(1, targetZRange(1) - processMargin), ...
    min(height, targetZRange(2) + processMargin) ...
];

z_start = processZRange(1);
z_end = processZRange(2);
n_sections = 5;
max_cap = 50;

% 横向只取中间有效区域估计 Power Doppler，避开左右两侧振铃和边缘伪影。
edgeMargin = 3;
x_start = min(edgeMargin + 1, width);
x_end = max(x_start, width - edgeMargin);
z_steps = round(linspace(z_start, z_end, n_sections + 1));

% 用候选点强度沿时间求平均能量，得到静态 Power Doppler 参考图。
PD_map = mean(MatInReduced.^2, 3);
PD_map = imgaussfilt(PD_map, 2);

roi_pd = PD_map(z_start:z_end, x_start:x_end);
pd_sorted = sort(roi_pd(roi_pd > 0), 'descend');

% 若 ROI 内存在有效能量点，则用排序后的较强能量作为补偿基准。
if ~isempty(pd_sorted)
    rank_ref = pd_sorted(min(10, length(pd_sorted)));
    factors = ones(n_sections, 1);
    z_pos = zeros(n_sections, 1);

    for s = 1:n_sections
        idx_r = z_steps(s):z_steps(s+1)-1;
        if isempty(idx_r)
            idx_r = z_steps(s);
        end

        sec_data = PD_map(idx_r, x_start:x_end);
        [m_val, m_idx] = max(sec_data(:));

        if m_val > 0
            [r_rel, ~] = ind2sub(size(sec_data), m_idx);
            z_pos(s) = idx_r(1) + r_rel - 1;
%             factors(s) = rank_ref / m_val;
            factors(s) = max(1, rank_ref / m_val);
        else
            z_pos(s) = mean(idx_r);
            factors(s) = 1;
        end
    end
    
    factors(end-3:end-1) = factors(end-3:end-1) * 2;
    factors = min(factors, max_cap);

    % 锁定最后一个控制点，保证增益曲线覆盖到缓冲区底部。
    z_pos(end) = z_end;
    all_z = [1; max(1, z_start - 1); z_pos; height];
    all_f = [1; 1; factors; factors(end)];

    [unique_z, u_idx] = unique(all_z);
    gain_curve_1d = interp1(unique_z, all_f(u_idx), 1:height, 'pchip')';
    gain_curve_1d = smoothdata(gain_curve_1d, 'gaussian', 5);
    gain_curve_1d = min(gain_curve_1d, max_cap);
    gain_curve_1d(~isfinite(gain_curve_1d)) = 1;

    fprintf('Gain debug | all_z / all_f: \n');
    disp([all_z(:), all_f(:)]);

    IntensityMatrix_Compensated = IntensityMatrix .* reshape(gain_curve_1d, [height, 1, 1]);
else
    IntensityMatrix_Compensated = IntensityMatrix;
end

% 增益估计完成后，只允许内部候选区域进入 Top-N 排序和定位。
candidateStart = candidateZRange(1);
candidateEnd = candidateZRange(2);
IntensityMatrix_Compensated(1:candidateStart-1, :, :) = 0;
IntensityMatrix_Compensated(candidateEnd+1:end, :, :) = 0;

% 继续屏蔽左右边缘列，减少 SVD 振铃造成的假点。
sideColumns = unique([1:min(edgeMargin, width), max(1, width-edgeMargin+1):width]);
IntensityMatrix_Compensated(:, sideColumns, :) = 0;

% 对补偿后的候选点逐帧排序，每帧只保留前 ULM.numberOfParticles 个候选点。
[tempMatrix, ~] = sort(reshape(IntensityMatrix_Compensated, [], numberOfFrames), 1, 'descend');
IntensityFinal = IntensityMatrix_Compensated - shiftdim(tempMatrix(ULM.numberOfParticles+1, :), -1);
MaskFinal = (IntensityFinal > 0) .* IntensityMatrix;
MaskFinal(~isfinite(MaskFinal)) = 0;

% 清理占内存较大的三维临时矩阵。
clear mask
clear IntensityMatrix IntensityMatrix_Compensated IntensityFinal
clear tempMatrix
clear PD_map roi_pd

% 将最终候选点从三维逻辑图转成坐标索引，供后续亚像素定位使用。
index_mask  = find(MaskFinal);
[index_mask_z, index_mask_x, index_numberOfFrames] = ind2sub([height, width, numberOfFrames], index_mask);
clear MaskFinal index_mask

%% 
averageXc = nan(1,size(index_mask_x,1),typename);% 创建一个全为 NaN的矩阵，所有候选微泡的 行索引
averageZc = nan(1,size(index_mask_z,1),typename);

for iscat=1:size(index_mask_z,1)
    % For each microbubble, create a 2D intensity matrix of the Region of interest defined by fwhm
    IntensityRoi = MatIn(index_mask_z(iscat)+vectfwhmz,index_mask_x(iscat)+vectfwhmx,index_numberOfFrames(iscat));%挖出一个 5*5的矩阵
%     IntensityRoi_norm = Mat_norm(index_mask_z(iscat)+vectfwhmz,index_mask_x(iscat)+vectfwhmx,index_numberOfFrames(iscat));%归一化数据Mat_norm
    % NLocal max
    % If there are too many localmax in the region of interest, the microbubble shape will be affected and the localization distorted.
    % In that case, we set averageZc, averageXc to NaN value.
    if nnz(imregionalmax(IntensityRoi))>ULM.parameters.NLocalMax
        continue
    end
    
        % 中心点至少要比周围背景平均值强约 9-2.8 / 6-2倍
    IntensityMax = MatIn(index_mask_z(iscat),index_mask_x(iscat),index_numberOfFrames(iscat));
    IntensityAverage = (sum(IntensityRoi(:))-IntensityMax)/(numel(IntensityRoi)-1);
    SNRLocal = 20*log10(IntensityMax / max(IntensityAverage, 1e-6));
%    
    if SNRLocal < 6
%         disp('deleted');
        continue;
    end
    
    % 判断当前候选亮点的局部形状，像不像一个理想的单个微泡高斯亮斑
    IntensityRoi_norm = (IntensityRoi - min(IntensityRoi(:)))./(max(IntensityRoi(:))-min(IntensityRoi(:)));
    if corr2(IntensityRoi_norm, Gauss_grid) < 0.65
        continue
    end
    
   
    switch ULM.LocMethod
        case 'curvefitting'
            [Zc,Xc,sigma] = curveFitting(IntensityRoi,vectfwhmz,vectfwhmx);
        case 'rs' 
            [Zc,Xc,sigma] = LocRadialSym(IntensityRoi,fwhmz,fwhmx);
    end
    
    % Store the final super-resolved position of the microbubble as its pixel position and an axial/lateral sub-pixel shift.
    averageZc(iscat) = Zc + index_mask_z(iscat);
    averageXc(iscat) = Xc + index_mask_x(iscat);% 最终位置 = 像素索引+亚像素偏移

    % Additional safeguards
    % sigma evaluates the size of the microbubble. If it appears to be too large, the microbubble can be removed (optional)
%     if or(sigma<0,sigma>25)
%         averageZc(iscat)=nan;
%          averageXc(iscat)=nan;
%         continue % 如果发现这种“畸形”的微泡，直接跳过，不保存它的坐标。
%     end

    % If the final axial/lateral shift is higher that the fwhmz,
    % localization has diverged and the microbubble is ignored.
    if or(abs(Zc)>fwhmz/2,abs(Xc)>fwhmx/2)% 偏移量超过了半径，不合格
        averageZc(iscat)=nan;
        averageXc(iscat)=nan;
%         M0(iscat) = nan;%矩
%         M2(iscat) = nan;
        continue
    end
end

keepIndex = ~isnan(averageXc);
ind = sub2ind([height,width,numberOfFrames],index_mask_z(keepIndex),index_mask_x(keepIndex),index_numberOfFrames(keepIndex));
clear index_mask_z index_mask_x IntensityRoi

%% BUILD MATTRACKING %%
% Creating the table which stores the high resolved microbubbles coordinates and the density value
MatTracking = zeros(nnz(keepIndex),4,typename);

MatTracking(:,1) = MatInReduced(ind);       % Initial intensity of the microbubble
MatTracking(:,2) = averageZc(keepIndex);    % Super-resolved axial coordinate
MatTracking(:,3) = averageXc(keepIndex);    % Super-resolved lateral coordinate
MatTracking(:,4) = index_numberOfFrames(keepIndex); % Frame number of the microbubble

MatTracking(:,2) = MatTracking(:,2) + delta(MatTracking(:,4),1);
MatTracking(:,3) = MatTracking(:,3) + delta(MatTracking(:,4),2);

% MatTracking(:,5) = M0(keepIndex);
% MatTracking(:,6) = M2(keepIndex);
clear averageXc averageZc index_numberOfFrames MatInReduced
end

function [Zc,Xc,sigma] = LocRadialSym(Iin,fwhm_z,fwhm_x)
%% function [Zc,Xc,sigma] = LocRadialSym(Iin,fwhm_z,fwhm_x)
    [Zc,Xc] = localizeRadialSymmetry(Iin,fwhm_z,fwhm_x);
    sigma = ComputeSigmaScat(Iin,Zc,Xc);
end


function sigma = ComputeSigmaScat(Iin,Zc,Xc)
%% This function will calculate the Gaussian width of the presupposed peak in the intensity, which we set as an estimate of the width of the microbubble
    [Nz,Nx] = size(Iin);
    Isub = Iin - mean(Iin(:));% 去除背景
    [pz,px] = meshgrid(1:Nx,1:Nz);% 生成像素网格
    zoffset = pz - Zc+(Nz)/2.0;%BH xoffset = px - xc;% 每一个像素点离微泡真实的“亚像素中心”到底有多远？
    xoffset = px - Xc+(Nx)/2.0;%BH yoffset = py - yc;
    r2 = zoffset.*zoffset + xoffset.*xoffset;% 距离的平方
    sigma = sqrt(sum(sum(Isub.*r2))/sum(Isub(:)))/2;  % second moment is 2*Gaussian width
end

function [Zc,Xc,sigma] = curveFitting(Iin,vectfwhm_z,vectfwhm_x)
%% function [Zc,Xc,sigma] = curveFitting(Iin,vectfwhm_z,vectfwhm_x)
% The ROI intensity is fitted with a theorical microbubble model.
    [meshX,meshZ] = meshgrid(vectfwhm_x,vectfwhm_z);
    meshIn = cat(3,meshX,meshZ);

    sigGauss_z = vectfwhm_z(end)*0+1;
    sigGauss_x = vectfwhm_x(end)*0+1;
    myGaussFunc = @(x_pos,mesh_pos)( exp(-(mesh_pos(:,:,1)-x_pos(1)).^2./(2*sigGauss_z^2) - (mesh_pos(:,:,2)-x_pos(2)).^2./(2*sigGauss_x^2)));
    OPTIONS = optimoptions('lsqcurvefit','StepTolerance',.01,'MaxIterations',5,'Display','off');

    % Gaussian fitting
    x_out = lsqcurvefit(myGaussFunc,[0 0],meshIn,double(Iin./max(Iin(:))),[],[],OPTIONS);
    Zc = x_out(2);
    Xc = x_out(1);
    sigma = ComputeSigmaScat(Iin,Zc,Xc);
end
