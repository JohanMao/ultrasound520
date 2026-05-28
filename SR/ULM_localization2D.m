function MatTracking = ULM_localization2D(MatIn,ULM,delta)
%ULM_LOCALIZATION2D 此处显示有关此函数的摘要
%   此处显示详细说明
fwhmz = ULM.fwhm(2);
fwhmx = ULM.fwhm(1);


vectfwhmz = -1*round(fwhmz/2):round(fwhmz/2);
vectfwhmx = -1*round(fwhmx/2):round(fwhmx/2);
[meshX,meshZ] = meshgrid(vectfwhmx,vectfwhmz);
sigGauss_z = vectfwhmz(end)*0+1;
sigGauss_x = vectfwhmx(end)*0+1;
Gauss_grid = exp(-(meshZ-0).^2./(2*sigGauss_z^2) - (meshX-0).^2./(2*sigGauss_x^2));
% [meshX,meshZ] = meshgrid(vectfwhmx,vectfwhmz);
% meshIn = cat(3,meshX,meshZ);
% 
% sigGauss_z = vectfwhm_z(end)*0+1;
% sigGauss_x = vectfwhm_x(end)*0+1;

[height,width,numberOfFrames]=size(MatIn);

MatIn = abs(MatIn);% Make sure you work with the intensity matrix
%在这里做复数转为实数，把变换放进去
Mat_norm = (MatIn-min(MatIn(:)))./(max(MatIn(:))-min(MatIn(:)));
info = whos('MatIn');typename = info.class;

if ~isfield(ULM,'LocMethod')
    ULM.LocMethod = 'curvefitting';
end

if ~isfield(ULM,'parameters')% Create an empty structure for parameters hosting
    ULM.parameters = struct();
end

if strcmp(ULM.LocMethod,'interp')
    if ~isfield(ULM.parameters,'InterpMethod')
        ULM.parameters.InterpMethod = 'spline';
    end
    if sum(strcmp(ULM.parameters.InterpMethod,{'bilinear','bicubic'}))
        warning('Faster but pixelated, Weighted Average will be faster and smoother.')
    end
end

if ~isfield(ULM.parameters,'NLocalMax')
    if fwhmz==3,ULM.parameters.NLocalMax = 2;
    else
        ULM.parameters.NLocalMax = 3;
    end
end

%% 
MatInReduced = zeros(height,width,numberOfFrames,typename);
MatInReduced(1+round(fwhmz/2)+1:height-round(fwhmz/2)-1,1+round(fwhmx/2)+1:width-round(fwhmx/2)-1,:) = MatIn(1+round(fwhmz/2)+1:height-round(fwhmz/2)-1, 1+round(fwhmx/2)+1:width-round(fwhmx/2)-1,:);
[height,width,numberOfFrames] = size(MatInReduced);

Mat2D = permute(MatInReduced, [1,3,2]); %so that all the frames are in columns
Mat2D = reshape(Mat2D,height*numberOfFrames,width);% Concatenate Matrix
mask2D = imregionalmax(Mat2D); clear Mat2D  % Perform imregionalmax
mask = reshape(mask2D,height,numberOfFrames,width);clear mask2D % reshape concatenated mask
mask = permute(mask,[1,3,2]); % so that we restore (z,x,t) table

IntensityMatrix = MatInReduced.*mask; %Values of intensities at regional maxima

% SELECTION OF MICROBUBBLES %%
% Only the first numberOfParticles highest local max will be kept for localization.
% Other local max will be considered as noise.
% Sort intensites in each frames, and store pixel coordinates
% At the end of this section, spatial and temporal coordinates microbubbles are
% stored into: index_mask_z, index_mask_x, index_numberOfFrames
[tempMatrix,~] = sort(reshape(IntensityMatrix,[],size(IntensityMatrix,3)),1,'descend');

% Remove the last kept intensity values to each frame. This means that you cannot fix an intensity threshold,
% we rely on number of particles. This is key for transparency/parallelization.
for i = 1:numberOfFrames
%    volumeBubble = IntensityMatrix(:,:,i);
%    volumeMax = max(volumeBubble(:));
%    threshIntensity = volumeMax * 10^(-20/20);
%    ULM.numberOfParticles = nnz(volumeBubble > threshIntensity);
% IntensityFinal(:,:,i) = IntensityMatrix(:,:,i) - ones(size(IntensityMatrix,1),size(IntensityMatrix,2)) .* reshape(tempMatrix( ULM.numberOfParticles+1,i),[1 1 numberOfFrames]);
    IntensityFinal(:,:,i) = IntensityMatrix(:,:,i) - ones(size(IntensityMatrix,1),size(IntensityMatrix,2)) .* tempMatrix(ULM.numberOfParticles+1,i);
end
clear tempMatrix
% Construction of the final mask with only the kept microbubbles low resolved and their associated intensity
MaskFinal = (mask.*IntensityFinal)>0;
MaskFinal(isnan(MaskFinal))=0;
MaskFinal = (MaskFinal>0).*IntensityMatrix;

% Preparing intensities and coordinates for further calculation of average, intensities etc...
index_mask  = find(MaskFinal);
[index_mask_z,index_mask_x,index_numberOfFrames]=ind2sub([height, width, numberOfFrames], index_mask);
clear mask IntensityFinal MaskFinal IntensityMatrix
clear index_mask
%% 
averageXc = nan(1,size(index_mask_z,1),typename);
averageZc = nan(1,size(index_mask_z,1),typename);

for iscat=1:size(index_mask_z,1)
    % For each microbubble, create a 2D intensity matrix of the Region of interest defined by fwhm
    IntensityRoi = MatIn(index_mask_z(iscat)+vectfwhmz,index_mask_x(iscat)+vectfwhmx,index_numberOfFrames(iscat));
    IntensityRoi_norm = Mat_norm(index_mask_z(iscat)+vectfwhmz,index_mask_x(iscat)+vectfwhmx,index_numberOfFrames(iscat));
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
    averageXc(iscat) = Xc + index_mask_x(iscat);

    % Additional safeguards
    % sigma evaluates the size of the microbubble. If it appears to be too large, the microbubble can be removed (optional)
    if or(sigma<0,sigma>50)
%         averageZc(iscat)=nan;
%         averageXc(iscat)=nan;
        continue
    end

    % If the final axial/lateral shift is higher that the fwhmz,
    % localization has diverged and the microbubble is ignored.
    if or(abs(Zc)>fwhmz/2,abs(Xc)>fwhmx/2)
        averageZc(iscat)=nan;
        averageXc(iscat)=nan;
%         M0(iscat) = nan;
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
    Isub = Iin - mean(Iin(:));
    [pz,px] = meshgrid(1:Nx,1:Nz);
    zoffset = pz - Zc+(Nz)/2.0;%BH xoffset = px - xc;
    xoffset = px - Xc+(Nx)/2.0;%BH yoffset = py - yc;
    r2 = zoffset.*zoffset + xoffset.*xoffset;
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

