clear all;
clc
% 
% datapath =  'H:\SR_data\12-Oct-2023\';
% datapath =  'E:\yj\lab\SR\PALA_data_InVivoMouseTumor\IQ1\';
% datapath =  'E:\yj\lab\motion correction\real_test\20241022-qy-rab-kid\temp\';
datapath =  'H:\data\20250411-skull\data\IQData\flow\';
format = '*.mat';
folder = dir([datapath,format]);
nameIQ = {folder.name};
Nbuffers = length(nameIQ);
load([datapath, nameIQ{2}]);
IQ = IQData;
% IQ = IregA;

% name = ['Rat1012_KSVD',num2str(k),'\'];
% name = ['Rat0427_KSVD_noS_singleThreshold\'];
% name = ['Tumor_KSVD_doubleThreshold\']; 
name = ['20250411flow\'];
localpath = ['H:\script\SR\Locals\',name];
trackspath = ['H:\script\SR\Tracks\',name];
savingpath = ['H:\script\SR\Mats\',name];
mkdir(localpath);
mkdir(trackspath);
mkdir(savingpath);
NFrames = size(IQ,3);
PData.Size = size(IQ);
PData.PDelta = [1 0 1]; 

UF.FrameRateUF = 1000;
UF.F0 = 15e6;
ULMparam = [200,10,1,8,50,249,4];

PData.Origin = [0 PData.Size(2)/2*PData.PDelta(2) 0];
framerate = UF.FrameRateUF;
%% ULM parameters
res = 10;

ULM = struct('numberOfParticles',200,...
    'res',res,...           
    'SVD_cutoff',round([10,NFrames]),...
    'max_linking_distance', 1,...
    'min_length',8,...
    'fwhm',[3 3],...
    'max_gap_closing',0,...
    'size',[PData.Size(1),PData.Size(2),NFrames],...
    'scale',[1 1 1/framerate],...
    'numberOfFrameProcessed',NFrames,...
    'interp_factor',1/res);

ULM.butter.CuttofFreq = [ULMparam(5) ULMparam(6)];
ULM.butter.samplingFreq = framerate;
[but_b, but_a] = butter(2,ULM.butter.CuttofFreq/(ULM.butter.samplingFreq/2),'bandpass');
ULM.parameters.NLocalMax = ULMparam(7);

pitch = 0.2/1000;
f0 = 10e6;
c = 1540;
lambda = c/f0;
Nelements = 128;
[M,N] = size(IQ(:,:,1));
% IQ = imresize(IQ,[M*2,N*2],'bilinear');
PData.Size = size(IQ);
Width = lambda * floor(N);
% Depth = lambda./2 * floor(M);   %为什么要除以2
Depth = lambda * floor(M);   %为什么要除以2
z = linspace(0,Depth,M).*1000;
x = linspace(0,Width,N).*1000;

% Nalgo = numel(listAlgo);

%% SVD filtering Noise
% ULM.SVD_cutoff = [2,20];
IQ_filt = SVDfilter(IQ,ULM.SVD_cutoff);
% IQ_filt = SVDfilter_temp(IQ);

% IQ_filt = filter(but_b, but_a,IQ_filt,[],3);
% IQ_filt(~isfinite(IQ_filt)) = 0;

% IQ_filt = IQ_filt./max(IQ_filt(:));
% IQ_filt = vesselness3D(IQ_filt,1:4,[1,1,1],0.75,true);

BullesSNR = abs(IQ_filt(:,:,1));
LocalMax = imregionalmax(BullesSNR);
ValMax = sort(BullesSNR(LocalMax),'descend');
noise = mean(BullesSNR(:));
thresh = noise * (10^(9/20));
mask = (ValMax >= thresh);

% ULM.numberOfParticles = round(mean(numberOfParticle(:)));
SNRmean = 20 * log10(mean(ValMax(1:10))/mean(BullesSNR(:)))
% clear BullesSNR LocalMax ValMax;

figure(1)
dB = 20 * log10(abs(IQ_filt));
dB = dB - max(dB(:));
imagesc(x,z,dB(:,:,1),[-40 0]),colormap gray;
colorbar;axis equal;axis tight;
xlabel('Axias Distance(mm)')
ylabel('Lateral Distance(mm)')
set(gca,'FontSize',15,'FontWeight','bold');
set(gca,'FontName','Times New Roman');
Bmode = dB(:,:,1);


figure(5)
PowDop = [];
for hhh=2:2
    path = [datapath, nameIQ{hhh}];
    tmp = load(path);
    IQ_filt = SVDfilter(tmp.IQData,ULM.SVD_cutoff);tmp = [];
    IQ_filt = filter(but_b,but_a,IQ_filt,[],3);
    IQ_filt(~isfinite(IQ_filt))=0;
    PowDop(:,:,end) = sqrt(sum(abs(IQ_filt).^2,3));

    pause(0.1)
end
im=imagesc(mean(PowDop,3).^(1/2));
axis image, colormap(gca,hot(128)),title(['Power Doppler'])
clbar = colorbar;
caxis([10 max(im.CData(:))*.9]);

xlabel('Axias Distance(mm)')
ylabel('Lateral Distance(mm)')
set(gca,'FontSize',15,'FontWeight','bold');
set(gca,'FontName','Times New Roman');

PD = mean(PowDop,3).^(1/2);
x1 = x;
z1 = z;

%% localize Data
clear Track_tot Track_tot_interp ProcessingTime IQ_filt IQ dB
fprintf('ULM PROCESSING\n');
t1 = tic;
hhh = 1;

% load('F:\SR_data_2\07-Mar-2023SonoVue_1SaveName')
% mcpath = [datapath,'mc\'];
mcpath = [datapath,'mc\'];

for hhh = 2:Nbuffers %1:Nbuffers
    fprintf('Processing bloc %d/%d\n',hhh,Nbuffers);

    path = [datapath, nameIQ{hhh}];
    tmp = load(path);
    
    mc = [mcpath,nameIQ{hhh}];
%     mctmp = load(mc);
    mctmp = {};
    mctmp.delta = zeros(800,2);
%     tmp.IQ = tmp.IQData;
    
    tmp.IQData(:,[1 end])=0;
    tmp.IQData([1 end],:)=0;
    IQ_filt = SVDfilter(tmp.IQData,ULM.SVD_cutoff);
%     IQ_filt = SVDfilter_s2d(tmp.IQData);
%     IQ_filt = SVDfilter_temp_st(tmp.IQ);
    
    IQ_filt = filter(but_b,but_a,IQ_filt,[],3);
    IQ_filt(~isfinite(IQ_filt))=0;
%     IQ_filt = BM3D_SR(IQ_filt);
    
    multiULM(IQ_filt,ULM,PData,mctmp.delta,'savingTrackfilename',[trackspath 'Tracks' num2str(hhh,'%.3d') '.mat'],...
    'savingLocalfilename',[localpath 'Locals' num2str(hhh,'%.3d') '.mat']);
end

% save([trackspath 'Tracks' num2str(1,'%.3d') '.mat'],'SNRmean','-append')
clear tmp IQ_filt

%%
MatOutSat = [];
MatOut= [];
MatOut_vel = MatOut;
MatOut_z = [];
MatOut_x = [];

fprintf('Building 2D ULM image\n');
index = 1;
MatOutBubble = [];
MatOutSat = [];

for hhh = 2:Nbuffers %1:Nbuffers
    trackLength = [];
    load([trackspath 'Tracks' num2str(hhh,'%.3d')],'Track_raw','Track_interp','ProTime');
    load([localpath 'Locals' num2str(hhh,'%.3d')]);
    aa = -PData(1).Origin([3 1])+[1 1]*1;
    bb = [1./PData(1).PDelta([3 1])];
    aa(3) = 0;
    bb(3:6) = 1;
    
    if size(cell2mat(Track_raw),2) == 4
        MatOutBubble(index,1) = 0;
        index = index + 1;
        continue;
    end
%     
%     [Track_raw,Track_interp] = kalman_filter(Track_raw,ULM); 
%     [Track_raw,Track_interp] = angle_constrain(Track_raw,Track_interp);
%     [Track_raw,Track_interp] = accela_constrain(Track_raw,Track_interp);
    
    Track_matout = Track_interp;
    Track_matout = cellfun(@(x) (x(:,[1 2 3 4 5 7]).*[bb]),Track_matout,'UniformOutput',0);
    Track_raw = cellfun(@(x) (x(:,[1 2 3 4 5 7]).*[bb]),Track_raw,'UniformOutput',0);
    for i = 1:length(Track_raw)
        singleTrack = Track_raw{i};
        trackLength(i) = length(singleTrack);
    end
    
    meanLength(hhh) = mean(trackLength);
    trackNumber(hhh) = length(Track_raw);
    trackBubble(hhh) = sum(trackLength(:))/size(MatTracking,1);
    
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
    
clear  Track_interp Track_count Track_matout 
save([savingpath 'MatOut_multi'],'MatOut','MatOut_vel','MatOutSat','ULM','PData','UF');

% trackLengthKSVD(k) = mean(meanLength);
% bubbleTrack(k) = sum(MatOut(:));
% trackedPercent(k) = mean(trackBubble);
% MatOutSatK{k} = MatOutSat;

% end

% load([savingpath 'MatOut_multi_slice0mm']);
%% 
% pitch = 0.3/1000;
% f0 = 7.5e6;
c = 1540;
lambda = c/f0;
Nelements = 128;
[M,N] = size(MatOut);
Width = pitch * floor(N/10);
Depth = lambda * floor(M/10);
z = linspace(0,Depth,M).*1000;
x = linspace(0,Width,N).*1000;
% MatZoom = [373 590 735 930];

figure(2)
IntPower = 1/2;SigmaGauss=0.3;
im=imagesc(x,z,imgaussfilt(MatOut.^IntPower,.01));axis image
if SigmaGauss>0,im.CData = imgaussfilt(im.CData,SigmaGauss);end

title('ULM intensity display')
colormap(gca,hot(128))
% clbar = colorbar;
caxis(caxis*.8)  % add saturation in image
% clbar.Label.String = 'number of counts';
% clbar.TickLabels = round(clbar.Ticks.^(1/IntPower),1);
% xlabel('\lambda');ylabel('\lambda') 
ca = gca;ca.Position = [.05 .05 .8 .9];

% im.CData = min(im.CData,10);caxis([0 10]) 
axis image
ylabel('Axias Distance(mm)')
xlabel('Lateral Distance(mm)')
% axis off
SR = imgaussfilt(MatOut.^IntPower,.01);


MatOut_vel = MatOut_vel./10;
figure(3)
vmax_disp  = ceil(quantile(MatOut_vel(abs(MatOut_vel)>0),.98)/10)*10;
IntPower = 1/2;
lambda = 1540/15e6 * 1e3;
ULM.SRscale = 10;

clf,set(gcf,'Position',[652 393 941 585]);
clbsize = [180,50];
Mvel_rgb = MatOut_vel/vmax_disp; % normalization
Mvel_rgb(1:clbsize(1),1:clbsize(2)) = repmat(linspace(1,0,clbsize(1))',1,clbsize(2)); % add velocity colorbar
Mvel_rgb = Mvel_rgb.^(1/2.5);Mvel_rgb(Mvel_rgb>1)=1;
Mvel_rgb = abs(Mvel_rgb);
Mvel_rgb = imgaussfilt(Mvel_rgb,.7);
Mvel_rgb = ind2rgb(round(Mvel_rgb*256),jet(256)); % convert ind into RGB

MatShadow = MatOut;MatShadow = MatShadow./max(MatShadow(:)*.3);MatShadow(MatShadow>1)=1;
MatShadow(1:clbsize(1),1:clbsize(2))=repmat(linspace(0,1,clbsize(2)),clbsize(1),1);
Mvel_rgb = Mvel_rgb.*(MatShadow.^IntPower);
Mvel_rgb = brighten(Mvel_rgb,.2);
BarWidth = round(1./(ULM.SRscale*lambda)); % 1 mm
Mvel_rgb(size(MatOut,1)-50+[0:3],60+[0:BarWidth],1:3)=1;
imshow(Mvel_rgb);axis on
title(['Velocity magnitude (0-' num2str(vmax_disp) 'mm/s)'])
ca = gca;ca.Position = [.05 .05 .8 .9];

figure(4)
MatOut_zdir = MatOut_vel;
velColormap = cat(1,flip(flip(hot(256),1),2),hot(256)); % custom velocity colormap
velColormap = velColormap(5:end-5,:); % remove white parts
IntPower = 1/2;
im=imagesc(x,z,(MatOut).^IntPower.*sign(imgaussfilt(MatOut_zdir,.3)));
im.CData = im.CData - sign(im.CData)/2;axis image
title(['ULM intensity display with axial flow direction'])
colormap(gca,velColormap)
caxis([-1 1]*max(caxis)*.8) % add saturation in image
clbar = colorbar;clbar.Label.String = 'Count intensity';
ca = gca;ca.Position = [.05 .05 .8 .9];


figure(6)
vmax_disp  = ceil(quantile(MatOut_int(abs(MatOut_int)>0),.98)/10)*10;
IntPower = 1/2;
lambda = 1;
ULM.SRscale = 10;
clf,set(gcf,'Position',[652 393 941 585]);
clbsize = [180,50];
Mvel_rgb = MatOut_int/vmax_disp; % normalization
Mvel_rgb(1:clbsize(1),1:clbsize(2)) = repmat(linspace(1,0,clbsize(1))',1,clbsize(2)); % add velocity colorbar
Mvel_rgb = Mvel_rgb.^(1/1.5);Mvel_rgb(Mvel_rgb>1)=1;
Mvel_rgb = abs(Mvel_rgb);
Mvel_rgb = imgaussfilt(Mvel_rgb,.5);
Mvel_rgb = ind2rgb(round(Mvel_rgb*256),hot(128)); % convert ind into RGB

MatShadow = MatOut;MatShadow = MatShadow./max(MatShadow(:)*.4);MatShadow(MatShadow>1)=1;
MatShadow(1:clbsize(1),1:clbsize(2))=repmat(linspace(0,1,clbsize(2)),clbsize(1),1);
Mvel_rgb = Mvel_rgb.*(MatShadow.^IntPower);
Mvel_rgb = brighten(Mvel_rgb,.3);
% BarWidth = round(1./(ULM.SRscale*lambda)); % 1 mm
% Mvel_rgb(size(MatOut,1)-50+[0:3],60+[0:BarWidth],1:3)=1;
imshow(Mvel_rgb);axis on
 
% 


% BarWidth = round(1./(ULM.scale(2)*lambda)); % 1 mm
% im.CData(size(im.CData,1)-2,3+[0:BarWidth])=max(caxis);


figure(7)
% t = linspace(1,round(length(MatOutSat)/20),length(MatOutBubble));
% plot(t, MatOutSat)
plot(smooth(MatOutSat,15));
xlabel('Data Batch');
ylabel('Number of Bubbles Tracked')


% 
% MatOut_structure = MatOut;
% MatOut_velocity = cat(3,MatOut_z,MatOut_x);
    