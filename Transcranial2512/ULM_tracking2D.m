function varargout = ULM_tracking2D( MatTracking,ULM,varargin )
%ULM_TRACKING2D 此处显示有关此函数的摘要
%   此处显示详细说明

%addpath('E:\yj\lab\SR\code\SVDfilter\SimpleTracker\SimpleTracker');

if nargin > 2 % Number of Arguments Input（输入参数的个数）
    mode = lower(varargin{1});%lower 变为小写
else
    mode = 'velocityinterp';
end

interp_factor = 1/ULM.max_linking_distance/ULM.res*.8;% 补帧因子
smooth_factor = 20;% 平滑因子
numberOfFrames = ULM.size(3);
FR = 1 / ULM.scale(3);

minFrame = min(MatTracking(:,4));
MatTracking(:,4) = MatTracking(:,4) - minFrame + 1;

index_Frame = arrayfun(@(i) find(MatTracking(:,4)==i),[1:numberOfFrames],'UniformOutput',false);% arrayfun批量处理 
% 得到了一个叫 index_Frame 的大柜子。拉开第 1 个抽屉（写成 index_Frame{1}），里面装的是 [5; 12; 108...]（意思是第 1 帧的微泡在大表里的第 5 行、12 行、108 行）。
Points = arrayfun(@(i) [MatTracking(index_Frame{i},2),MatTracking(index_Frame{i},3)],...
    [1:numberOfFrames],'UniformOutput',false);% 拿“真实坐标”
% ,MatTracking(index_Frame{i},5),MatTracking(index_Frame{i},6)
debug = false;
[ Simple_Tracks,Adjacency_Tracks ] = simpletracker(Points,...
    'MaxLinkingDistance', ULM.max_linking_distance, ...
    'MaxGapClosing', ULM.max_gap_closing, ...
    'Debug', debug);
n_tracks = numel(Simple_Tracks);% 共找出了多少条微泡轨迹
all_points = vertcat(Points{:});% 垂直拼接（Vertical Concatenate）

count=1;
Tracks_raw = {};

for i_track = 1:n_tracks
    track_id = Adjacency_Tracks{i_track};
    idFrame = MatTracking(track_id,4);
    intBubble = MatTracking(track_id,1);
    track_points = cat(2,all_points(track_id,:),idFrame,intBubble);
    
    if length(track_points(:,1)) > ULM.min_length
        
%         % 加速度约束 (Acceleration Constraint)
%         v_vec = diff(track_points(:,1:2), 1, 1); 
%         v = sqrt(sum(v_vec.^2, 2));
%         v_trim_mean = trimmean(v, 20);
% 
%         a_thr = 1.5 * max(v_trim_mean, eps) * FR;
%         a = abs(diff(v)) * FR;
% 
%         if sum(a > a_thr) > max(2, round(0.1 * numel(a)))
%             continue;
%         end
        
%         % 方向约束 (Direction Constraint), 微泡不能发生超过 90 度的锐角折返 
%         is_valid_direction = true;
%         for j = 1:(size(v_vec, 1) - 1)
%             v1 = v_vec(j, :);
%             v2 = v_vec(j+1, :);
%             cos_theta = dot(v1, v2) / (norm(v1) * norm(v2) + eps);
%             if cos_theta < 0 % 余弦值小于0意味着夹角大于90度
%                 is_valid_direction = false;
%                 break;
%             end
%         end
%         if ~is_valid_direction
%             continue; 
%         end
        
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


