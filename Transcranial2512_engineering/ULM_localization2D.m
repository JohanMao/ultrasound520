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
% % --- 调试：查看 5-10 行是否有局部极大值点 ---
% % 1. 截取第 5 到 10 行的所有数据
% mask_slice = mask(50:100, :, :); 
% 
% % 2. 找到所有值为 1 (极大值) 的线性索引
% idx = find(mask_slice); 
% 
% if isempty(idx)
%     disp('--- 警告：第 5-10 行没有发现任何局部极大值点！---');
% else
%     % 3. 将线性索引转回 (z, x, frame) 坐标
%     [z_rel, x, f] = ind2sub(size(mask_slice), idx);
%     
%     % 4. 把相对坐标 z_rel 转回原始行号 (5-10)
%     z_abs =z_rel + (50 - 1);
%     
%     % 5. 组合成坐标表 [行号, 列号, 帧号]
%     found_points = [z_abs, x, f];
%     
%     fprintf('--- 在第 5-10 行共抓到 %d 个候选点 ---\n', size(found_points, 1));
%     % 显示前 10 个点看看情况
%     disp('前 10 个候选点坐标 [Row, Col, Frame]:');
%     disp(found_points(1:min(10, end), :));
% end

IntensityMatrix = MatInReduced.*mask; %Values of intensities at regional maxima保留那些被标记为候选点的像素亮度值

% %% 增益补偿模块：仅针对第 5 行至第 27 行
% n_sections = 10; 
% z_start_idx = 47; 
% z_end_idx = 107;
% range_len = z_end_idx - z_start_idx + 1;
% 
% % 
% z_steps = round(linspace(z_start_idx, z_end_idx, n_sections + 1));
% gain_factors = ones(n_sections, numberOfFrames);
% gain_curve = ones(height, numberOfFrames); % 默认全1，保证区间外不增强
% 
% for i = 1:numberOfFrames
%     frame_data = IntensityMatrix(:,:,i);
%     global_max = max(frame_data(:));
%     if global_max == 0, continue; end 
% %     all_pos = frame_data(frame_data > 0); % 仅提取该帧中有亮度的候选点
% %     if isempty(all_pos), continue; end
% %     
% %     % 取全图前 1% 亮点的平均值（比 max 稳健得多）
% %     sorted_all = sort(all_pos, 'descend');
% %     n_all = max(1, round(0.001 * length(sorted_all))); 
% %     global_ref = mean(sorted_all(1:n_all));
%     
%     
%     for s = 1:n_sections
%         section_idx = z_steps(s):z_steps(s+1) - 1;
%         section_max = max(max(frame_data(section_idx, :)));
% %         section_mean = mean(frame_data(section_idx, :), 'all');
%         
%         if section_max > 0
%             % 无上限补偿
%             gain_factors(s, i) = global_max / section_max;
%         end
%     end
%     
%     z_centers = (z_steps(1:end-1) + z_steps(2:end)) / 2;
%     care_range = z_start_idx:z_end_idx;
%     gain_curve(care_range, i) =  0.5 * interp1(z_centers, gain_factors(:, i), care_range, 'pchip', 1);% 三次Hermite插值
%     
% end
% 
% % 应用补偿
% IntensityMatrix_Compensated = IntensityMatrix .* reshape(gain_curve, [height, 1, numberOfFrames]);
% 
% % IntensityMatrix_Compensated(50:70, :, :) = 0;
% [tempMatrix, ~] = sort(reshape(IntensityMatrix_Compensated, [], numberOfFrames), 1, 'descend');

% %% 改进的增益补偿模块：基于 Max 且分两大段处理
% % 第一部分：47-91 行 
% n_sections_1 = 5;
% z_start_1 = 47; 
% z_end_1   = 91;
% 
% % 第二部分：92-107 行 (分成 4 段)
% n_sections_2 = 4;
% z_start_2 = 92;
% z_end_2   = 107;
% 
% z_steps_1 = round(linspace(z_start_1, z_end_1, n_sections_1 + 1));
% z_steps_2 = round(linspace(z_start_2, z_end_2, n_sections_2 + 1));
% 
% gain_curve = ones(height, numberOfFrames); 
% 
% for i = 1:numberOfFrames
%     frame_data = IntensityMatrix(:,:,i);
%     % 全局最大值基准
%     global_max = max(frame_data(:));
%     if global_max == 0, continue; end
%     all_pos = frame_data(frame_data > 0);
%     if isempty(all_pos), continue; end
%     % --- 核心改进：提取全图排序第 200 名的强度作为基准 ---
%     s_all = sort(all_pos, 'descend');
%     rank_ref1 = s_all(min(150, length(s_all))); 
%     rank_ref2 = s_all(min(100, length(s_all)));
% %     if rank_ref == 0, rank_ref = max(s_all); end
% 
%     % --- 第一部分计算：47-91 行 ---
%     factors_1 = ones(n_sections_1, 1);
%     for s = 1:n_sections_1
%         idx = z_steps_1(s):z_steps_1(s+1)-1;
%         section_max = max(max(frame_data(idx, :)));
%         if section_max > 0
%             % 无上限补偿：直接拉齐到全局最大值
%             factors_1(s) = rank_ref1 / section_max;
%         end
%     end
%     
%     % --- 第二部分计算：92-107 行 (分成 4 段) ---
%     factors_2 = ones(n_sections_2, 1);
%     for s = 1:n_sections_2
%         idx = z_steps_2(s):z_steps_2(s+1)-1;
%         section_max = max(max(frame_data(idx, :)));
%         if section_max > 0
%             % 无上限补偿
%             factors_2(s) = rank_ref2 / section_max;
%         end
%     end
% 
%     % --- 分段插值：使用 'extrap' 确保覆盖到 107 行 ---
%     % 第一段插值 (47-91)
%     z_c1 = (z_steps_1(1:end-1) + z_steps_1(2:end)) / 2;
%     range_1 = z_start_1:z_end_1; %动态控制点定位
%     gain_curve(range_1, i) = interp1(z_c1, factors_1, range_1, 'pchip', 'extrap');
%     
%     % 第二段插值 (92-107)
%     z_c2 = (z_steps_2(1:end-1) + z_steps_2(2:end)) / 2;
%     range_2 = z_start_2:z_end_2;
%     gain_curve(range_2, i) = interp1(z_c2, factors_2, range_2, 'linear', 'extrap');
%     
%     % 可选：平滑两段之间的接缝
%     gain_curve(47:107, i) = smoothdata(gain_curve(47:107, i), 'gaussian', 3);
% end
% 
% IntensityMatrix_Compensated = IntensityMatrix .* reshape(gain_curve, [height, 1, numberOfFrames]);
% % 应用补偿
% [tempMatrix, ~] = sort(reshape(IntensityMatrix_Compensated, [], numberOfFrames), 1, 'descend');
%%
%% 改进的增益补偿模块：基于 Max 且分两大段处理
% 第一部分：47-91 行 
n_sections_1 = 5;
z_start_1 = 70; 
z_end_1   = 85;

% 第二部分：92-107 行 (分成 4 段)
n_sections_2 = 7;
z_start_2 = 86;
z_end_2   = 107;

z_steps_1 = round(linspace(z_start_1, z_end_1, n_sections_1 + 1));
z_steps_2 = round(linspace(z_start_2, z_end_2, n_sections_2 + 1));
gain_curve = ones(height, numberOfFrames); 

% 改进的动态控制点增益补偿

for i = 1:numberOfFrames
    frame_data = IntensityMatrix(:,:,i);
    global_max = max(frame_data(:));
    if global_max == 0, continue; end
    all_pos = frame_data(frame_data > 0);
    if isempty(all_pos), continue; end
    % --- 核心改进：提取全图排序第 200 名的强度作为基准 ---
    s_all = sort(all_pos, 'descend');
%     rank_ref = s_all(min(70, length(s_all)));
    rank_ref1 = s_all(min(75, length(s_all))); 
    rank_ref2 = s_all(min(50, length(s_all)));
    % --- 处理第一段 (47-91) ---
    factors_1 = ones(n_sections_1, 1);
    z_pos_1 = zeros(n_sections_1, 1);
    for s = 1:n_sections_1
        idx_r = z_steps_1(s):z_steps_1(s+1)-1;
        sec_data = frame_data(idx_r, :);
        [m_val, m_idx] = max(sec_data(:));% 数值、线性索引
        if m_val > 0
            [r_rel, ~] = ind2sub(size(sec_data), m_idx);
            z_pos_1(s) = idx_r(1) + r_rel - 1;
            factors_1(s) = rank_ref1 / m_val;
        else
            z_pos_1(s) = (idx_r(1) + idx_r(end))/2;
            factors_1(s) = 1;
        end
    end
    % 关键：确保控制点不重叠且有序（防止极端情况）
%     z_pos_1 = sort(unique(z_pos_1)); 
%     if length(z_pos_1) < n_sections_1, z_pos_1 = linspace(z_start_1, z_end_1, n_sections_1); end % 兜底
    gain_curve(z_start_1:z_end_1, i) = interp1(z_pos_1, factors_1, z_start_1:z_end_1, 'pchip', 'extrap');

    % --- 处理第二段 (92-107) ---
    factors_2 = ones(n_sections_2, 1);
    z_pos_2 = zeros(n_sections_2, 1);
    for s = 1:n_sections_2
        idx_r = z_steps_2(s):z_steps_2(s+1)-1;
        sec_data = frame_data(idx_r, :);
        [m_val, m_idx] = max(sec_data(:));
        if m_val > 0
            [r_rel, ~] = ind2sub(size(sec_data), m_idx);
            z_pos_2(s) = idx_r(1) + r_rel - 1;
            factors_2(s) = rank_ref2 / m_val;
        else
            z_pos_2(s) = (idx_r(1) + idx_r(end))/2;
            factors_2(s) = 1;
        end
    end
%     z_pos_2 = sort(unique(z_pos_2));
    % 强制最后一个点在 107 行，解决你之前的边缘断层
    z_pos_2(end) = z_end_2; 
    gain_curve(z_start_2:z_end_2, i) = interp1(z_pos_2, factors_2, z_start_2:z_end_2, 'pchip', 'extrap');
    gain_curve(47:107, i) = smoothdata(gain_curve(47:107, i), 'gaussian', 3);
end

IntensityMatrix_Compensated = IntensityMatrix .* reshape(gain_curve, [height, 1, numberOfFrames]);
[tempMatrix, ~] = sort(reshape(IntensityMatrix_Compensated, [], numberOfFrames), 1, 'descend');

%%

% SELECTION OF MICROBUBBLES %%
% Only the first numberOfParticles highest local max will be kept for localization.
% Other local max will be considered as noise.
% Sort intensites in each frames, and store pixel coordinates
% At the end of this section, spatial and temporal coordinates microbubbles are
% stored into: index_mask_z, index_mask_x, index_numberOfFrames

% [tempMatrix,~] = sort(reshape(IntensityMatrix,[],size(IntensityMatrix,3)),1,'descend');% 从高到低进行排名，1代表沿着第一维，每一帧里全场最亮的那个点的亮度，~ 丢弃了索引，tempMatrix是具体的亮度值
%它的第 1 行就是每一帧里全场最亮的那个点的亮度；第 2 行是每一帧次亮的亮度。
% Remove the last kept intensity values to each frame. This means that you cannot fix an intensity threshold,
% we rely on number of particles. This is key for transparency/parallelization.
for i = 1:numberOfFrames
%    volumeBubble = IntensityMatrix(:,:,i);
%    volumeMax = max(volumeBubble(:));
%    threshIntensity = volumeMax * 10^(-20/20);
%    ULM.numberOfParticles = nnz(volumeBubble > threshIntensity);
% IntensityFinal(:,:,i) = IntensityMatrix(:,:,i) - ones(size(IntensityMatrix,1),size(IntensityMatrix,2)) .* reshape(tempMatrix( ULM.numberOfParticles+1,i),[1 1 numberOfFrames]);
%    IntensityFinal(:,:,i) = IntensityMatrix(:,:,i) - ones(size(IntensityMatrix,1),size(IntensityMatrix,2)) .* tempMatrix(ULM.numberOfParticles+1,i);
%找出前numberOfFrames个点的强度
   IntensityFinal(:,:,i) = IntensityMatrix_Compensated(:,:,i) - tempMatrix(ULM.numberOfParticles+1, i);
end
%       IntensityFinal(:,:,i) = IntensityMatrix(:,:,i) - ones(size(IntensityMatrix,1),size(IntensityMatrix,2)) .* tempMatrix(ULM.numberOfParticles+1,i);% 找出前numberOfParticles个点的强度
% end
% clear tempMatrix
% Construction of the final mask with only the kept microbubbles low resolved and their associated intensity
% MaskFinal = (mask.*IntensityFinal)>0;
% MaskFinal(isnan(MaskFinal))=0;% 找矩阵中所有的NaN并将其设为 0。
% MaskFinal = (MaskFinal>0).*IntensityMatrix;
% MaskFinal = (MaskFinal>0).*IntensityMatrix_Compensated;
MaskFinal = (IntensityFinal > 0) .* IntensityMatrix_Compensated;
% MaskFinal = (IntensityFinal > 0) .* IntensityMatrix;
MaskFinal(~isfinite(MaskFinal)) = 0;

% Preparing intensities and coordinates for further calculation of average, intensities etc...
index_mask  = find(MaskFinal);% 把三维的数据块拉成一条直线，告诉你哪几个位置上有气泡
[index_mask_z,index_mask_x,index_numberOfFrames]=ind2sub([height, width, numberOfFrames], index_mask);% 把上面那个线性索引转换回
clear mask IntensityFinal MaskFinal IntensityMatrix_Compensated
clear index_mask
%% 
averageXc = nan(1,size(index_mask_x,1),typename);% 创建一个全为 NaN的矩阵，所有候选微泡的 行索引
averageZc = nan(1,size(index_mask_z,1),typename);

for iscat=1:size(index_mask_z,1)
    % For each microbubble, create a 2D intensity matrix of the Region of interest defined by fwhm
    IntensityRoi = MatIn(index_mask_z(iscat)+vectfwhmz,index_mask_x(iscat)+vectfwhmx,index_numberOfFrames(iscat));%挖出一个 5*5的矩阵
    IntensityRoi_norm = Mat_norm(index_mask_z(iscat)+vectfwhmz,index_mask_x(iscat)+vectfwhmx,index_numberOfFrames(iscat));%归一化数据Mat_norm
    % NLocal max
    % If there are too many localmax in the region of interest, the microbubble shape will be affected and the localization distorted.
    % In that case, we set averageZc, averageXc to NaN value.
%     if nnz(imregionalmax(IntensityRoi))>ULM.parameters.NLocalMax
%         continue
%     end
    
%     IntensityRoi_norm = (IntensityRoi - min(IntensityRoi(:)))./(max(IntensityRoi(:))-min(IntensityRoi(:)));
% %     
%     if corr2(IntensityRoi_norm, Gauss_grid) < 0.65
%         continue
%     end
% % % %     
%    IntensityMax = MatIn(index_mask_z(iscat),index_mask_x(iscat),index_numberOfFrames(iscat));
%    IntensityAverage = (sum(IntensityRoi(:))-IntensityMax)/(numel(IntensityRoi)-1);
%    SNRLocal = 20*log10(IntensityMax/IntensityAverage);
% %    
%    if SNRLocal < 9
% %        disp('deleted');
%        continue;
%    end
   
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
%     if or(sigma<0,sigma>10)
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