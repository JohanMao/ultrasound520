function [Track_raw,Track_interp,varargout] = multiULM(IQ,ULM,PData,delta,varargin)

if ~isfield(ULM, 'parameters')% 检查结构体 ULM 中是否包含名为 'parameters' 的字段
    ULM.parameters = struct();
end

if ~isfield(ULM.parameters, 'NLocalMax')
    ULM.parameters.NLocalMax = 4;
end

tmp = strcmpi(varargin, 'tracking');
if any(tmp),tracking = varargin{find(tmp) + 1};else,tracking=1;end

tmp = strcmpi(varargin, 'savingTrackfilename');
if any(tmp),savingTrackingFileName = varargin{find(tmp)+1};SaveTrackData = true;else,SaveTrackData = false;end

tmp = strcmpi(varargin, 'savingLocalfilename');
if any(tmp),savingLocalFileName = varargin{find(tmp)+1};SaveLocalData = true;else,SaveLocalData = false;end

ULM.max_linking_distance = ULM.max_linking_distance * PData.PDelta(3); % 从像素单位转化为波长
ULM.LocMethod = 'rs';  %径向对称法（Radial Symmetry）

localizationTimer = tic;
[MatTracking] = ULM_localization2D(IQ,ULM,delta);
localizationTime = toc(localizationTimer);

if SaveLocalData
    save(savingLocalFileName,'MatTracking')% 把MatTracking存起来
end

% ProcessingTime = 0;
% load(savingLocalFileName);
% MatTracking(:,2) = MatTracking(:,2) + delta(MatTracking(:,4),1);
% MatTracking(:,3) = MatTracking(:,3) + delta(MatTracking(:,4),2);

%%
Track_raw = {};
MatTracking(:,2:3) = (MatTracking(:,2:3) - repmat([1 1],[size(MatTracking,1),1]))...
    .* repmat(PData.PDelta([3 1]),[size(MatTracking,1),1]) + repmat([PData.Origin(3) PData.Origin(1)],[size(MatTracking,1),1]);
% MatTracking(:,2:3) = (MatTracking(:,2:3) - 1) .* PData.PDelta([3 1]) + PData.Origin([3 1]);% 物理坐标 = （像素网格位置-1）*单个网格的真实物理大小+ 真实世界原点
if tracking
    trackingTimer = tic;
    [Track_interp,Track_raw] = ULM_tracking2D(double(MatTracking),ULM,'velocityinterp');
    trackingTime = toc(trackingTimer);
    Track_interp = cellfun(@single,Track_interp,'UniformOutput',false);
    Track_raw = cellfun(@single,Track_raw,'UniformOutput',false);%把所有轨迹数据从‘高精度、大内存占用’模式转换成‘够用精度、小内存占用’模式。
    ProcessingTime = localizationTime + trackingTime;
    fprintf('Time | localization: %.3f s, tracking: %.3f s\n', localizationTime, trackingTime);
else
    Track_raw = single(MatTracking);% 把所有的 double 重新转回 single
    Track_interp = [];
    ProcessingTime = localizationTime;
    fprintf('Time | localization: %.3f s\n', localizationTime);
end

    
if nargout ==3
    varargout{1} = ProcessingTime;% 输出参数计数器 (Number of Arguments Out)
end

%% Save data in .mat file
if SaveTrackData
    fprintf('saving... ')
    ProTime = ProcessingTime;
    save(savingTrackingFileName,'Track_interp','Track_raw','ProTime','ULM','PData');
%     'Track_raw'
end
fprintf('end.\n')
end
