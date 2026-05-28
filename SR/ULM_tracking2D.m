function varargout = ULM_tracking2D( MatTracking,ULM,varargin )
%ULM_TRACKING2D 此处显示有关此函数的摘要
%   此处显示详细说明

addpath('E:\yj\lab\SR\code\SVDfilter\SimpleTracker\SimpleTracker');

if nargin > 2
    mode = lower(varargin{1});
else
    mode = 'velocityinterp';
end

interp_factor = 1/ULM.max_linking_distance/ULM.res*.8;

smooth_factor = 20;
numberOfFrames = ULM.size(3);

minFrame = min(MatTracking(:,4));
MatTracking(:,4) = MatTracking(:,4) - minFrame + 1;
index_Frame = arrayfun(@(i) find(MatTracking(:,4)==i),[1:numberOfFrames],'UniformOutput',false);
Points = arrayfun(@(i) [MatTracking(index_Frame{i},2),MatTracking(index_Frame{i},3)],...
    [1:numberOfFrames],'UniformOutput',false);
% ,MatTracking(index_Frame{i},5),MatTracking(index_Frame{i},6)
debug = false;
tic
[ Simple_Tracks,Adjacency_Tracks ] = simpletracker(Points,...
    'MaxLinkingDistance', ULM.max_linking_distance, ...
    'MaxGapClosing', ULM.max_gap_closing, ...
    'Debug', debug);
toc
n_tracks=numel(Simple_Tracks);
all_points = vertcat(Points{:});

count=1;
Tracks_raw = {};

for i_track = 1:n_tracks
    track_id = Adjacency_Tracks{i_track};
    idFrame = MatTracking(track_id,4);
    intBubble = MatTracking(track_id,1);
    track_points = cat(2,all_points(track_id,:),idFrame,intBubble);
    if length(track_points(:,1))>ULM.min_length
        Tracks_raw{count}=track_points;
        count=count+1;
    end
end

if count==1
    disp(['Was not able to find tracks at ',num2str(minFrame)]);
    Tracks_out{1}=[0,0,0,0];varargout{1}=Tracks_out;
    if nargout>1,varargout{2} = Tracks_out;end
    return
end

Tracks_out = {};

for i_track = 1:size(Tracks_raw,2)
    track_points=double(Tracks_raw{1,i_track});
    xi=track_points(:,2);
    zi=track_points(:,1);
    ii = track_points(:,4);
    
    TimeAbs=(0:(length(zi)-1))*ULM.scale(3);
    frame = track_points(:,3);
    % Interpolation of spatial and time components
    zu=interp1(1:length(zi),smooth(zi),1:interp_factor:length(zi));
    xu=interp1(1:length(xi),smooth(xi),1:interp_factor:length(xi));
    iiu=interp1(1:length(ii),smooth(ii),1:interp_factor:length(ii));
    
    TimeAbs_interp = interp1(1:length(TimeAbs),TimeAbs,1:interp_factor:length(TimeAbs));
    frame_interp = interp1(1:length(frame),frame,1:interp_factor:length(frame));
    % Velocity
    vz = diff(zi')./diff(TimeAbs);vz=[vz(1),vz];
    vx = diff(xi')./diff(TimeAbs);vx=[vx(1),vx];
    vzu=diff(zu)./diff(TimeAbs_interp);vzu=[vzu(1),vzu];
    vxu=diff(xu)./diff(TimeAbs_interp);vxu=[vxu(1),vxu];
    
    if length(zi)>ULM.min_length
        Tracks_raw_v{i_track,1}=single(cat(2,zi,xi,vz',vx',TimeAbs',frame,ii));
        Tracks_out{i_track,1}=single(cat(2,zu',xu',vzu',vxu',TimeAbs_interp',frame_interp',iiu')); %position / velocity / timeline
    end
end

Tracks_out = Tracks_out(~cellfun('isempty',Tracks_out));
varargout{1}=Tracks_out;
varargout{2}=Tracks_raw_v;


