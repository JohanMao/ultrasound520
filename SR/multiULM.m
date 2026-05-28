function [Track_raw,Track_interp,varargout] = multiULM(IQ,ULM,PData,delta,varargin)

if ~isfield(ULM, 'parameters')
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

ULM.max_linking_distance = ULM.max_linking_distance * PData.PDelta(3);

ULM.LocMethod = 'rs';
t0 = tic;
[MatTracking] = ULM_localization2D(IQ,ULM,delta);
ProcessingTime = toc(t0);

if SaveLocalData
    save(savingLocalFileName,'MatTracking')
end

% ProcessingTime = 0;
% load(savingLocalFileName);
% MatTracking(:,2) = MatTracking(:,2) + delta(MatTracking(:,4),1);
% MatTracking(:,3) = MatTracking(:,3) + delta(MatTracking(:,4),2);

%%
Track_raw = {};
MatTracking(:,2:3) = (MatTracking(:,2:3) - repmat([1 1],[size(MatTracking,1),1]))...
    .* repmat(PData.PDelta([3 1]),[size(MatTracking,1),1]) + repmat([PData.Origin(3) PData.Origin(1)],[size(MatTracking,1),1]);
if tracking
        [Track_interp,Track_raw] = ULM_tracking2D(double(MatTracking),ULM,'velocityinterp');
    else
         Track_raw = single(MatTracking);
        Track_interp = [];
end

    Track_interp = cellfun(@single,Track_interp,'UniformOutput',false);
    Track_raw = cellfun(@single,Track_raw,'UniformOutput',false);

    
if nargout ==3
    varargout{1} = ProcessingTime;
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
