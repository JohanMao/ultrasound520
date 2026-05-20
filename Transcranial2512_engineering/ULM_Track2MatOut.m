function [MatOut,varargout] = ULM_Track2MatOut(Tracks,sizeOut,varargin)
%% function MatOut = ULM_Track2MatOut(Tracks,sizeOut,varargin)
% This function returns a density image from all tracks. Each track is projected on a
% grid image. Pixel intensity represents the number of tracks crossing this pixel.
% INPUTS:
%       - Tracks is the list of tracks. Positions in x and y are indexes of pixels in the rendering grid (ie. 1/res). The conversion
%       has to be done before sending tracks to this function.
%       - sizeOut is the size of the grid where tracks are projected.
%       - 'mode' (in varargin):
%               - 2D_velnorm : mean velocity image
%               - 2D_velmean : mean velocity per track
%               - 2D_vel_z : to reconstruct velocity postive upward, and negative if downward
%           if empty, 'mode' will be set to 2D_tracks or 2D_allin based on whether input 'Tracks' are cells or not.
% OUTPUTS:
%       - MatOut is the density map of size sizeOut (generally 10x10 times the size of input IQ)
%
% Created by Arthur Chavignon, 11/04/19
%
% DATE 2020.07.22 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code published in 2020 article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.
% Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018

%% Input positions are in pixels units (for example z=1.2,x=3.54), then converted in index
% in the grid by rounding value (z=1,x=4). The pixel of position (1,4) will be increment by 1.

%% select mode
tmp = strcmpi(varargin,'mode');
if any(tmp),mode= varargin{find(tmp)+1};else, mode=[];end

if isempty(mode)
    if iscell(Tracks)
        mode = '2D_tracks';
    else
        mode = '2D_allin';
    end
end
res = varargin{4}.res;
%% Build MatOut
MatOut = zeros(sizeOut); % Initialization of the density image
if ismember(mode,{'2D_vel_z','2D_velnorm','2D_velmean'})
    MatOut_vel = zeros(sizeOut);
    MatOut_vel_z = zeros(sizeOut);
    MatOut_vel_x = zeros(sizeOut);
    MatOut_int = zeros(sizeOut);
end

switch mode
    case '2D_allin'
        %% 2D case - 1 cell
%         disp(['Building 2D ULM image (size : ' num2str(sizeOut) ').'])
        % round pixel position into [z,x] index
        pos_z_round = round(Tracks(:,1));
        pos_x_round = round(Tracks(:,2));

        % remove out of grid bubbles (ie. the grid is too small)
        outP = zeros(numel(pos_z_round),1);
        outP = or(outP,pos_z_round<1);outP = or(outP,pos_z_round>sizeOut(1));
        outP = or(outP,pos_x_round<1);outP = or(outP,pos_x_round>sizeOut(2));

        ind=(sub2ind(sizeOut,pos_z_round(~outP),pos_x_round(~outP)));

        for j=1:size(ind,1)
            % increment the intensity count of the pixel crossed by a track by +1
            MatOut(ind(j))=MatOut(ind(j))+1;
        end

    case '2D_tracks'
        %% 2D case - 1 cell per tracks
%         disp(['Building 2D ULM image (size : ' num2str(sizeOut) ').'])
        for itrack = 1:numel(Tracks)
            % round pixel position into [z,x] index
            pos_z_round = round(Tracks{itrack}(:,1));
            pos_x_round = round(Tracks{itrack}(:,2));

            % remove out of grid bubbles (ie. the grid is too small)
            outP = zeros(numel(pos_z_round),1);
            outP = or(outP,pos_z_round<1);outP = or(outP,pos_z_round>sizeOut(1));
            outP = or(outP,pos_x_round<1);outP = or(outP,pos_x_round>sizeOut(2));

            ind=(sub2ind(sizeOut,pos_z_round(~outP),pos_x_round(~outP)));
            % Tracks are counted only once per pixel, the unique keeps only 1 point per pixel
            ind = unique(ind);

            for j=1:size(ind,1)
                % increment the intensity count of the pixel crossed by a track by +1
                MatOut(ind(j))=MatOut(ind(j))+1;
            end
        end

    case {'2D_vel_z','2D_velnorm','2D_velmean'}
        %% 2D case - 1 cell per tracks
        % Creates mean velocity matrix
%         disp(['Building 2D ULM image (size : ' num2str(sizeOut) ').'])
        for itrack = 1:numel(Tracks)
            % round pixel position into [z,x] index
            pos_z_round = round(Tracks{itrack}(:,1).*res);
            pos_x_round = round(Tracks{itrack}(:,2).*res);
            
            % intBubble = smooth(Tracks{itrack}(:,6),10);
            % 修改后的代码 (使用内置的移动平均平滑)
            intBubble = smoothdata(Tracks{itrack}(:,6), 'movmean', 10);
            if strcmp(mode,'2D_velmean')
                velnorm = Tracks{itrack}(:,3);
            else
%                  velnorm = smooth(vecnorm(Tracks{itrack}(:,3:4),2,2),10); % velocity is smoothed
                 velnorm = smoothdata(vecnorm(Tracks{itrack}(:,3:4),2,2), 'movmean', 10);
%                  velnorm_z = smooth(Tracks{itrack}(:,3),10);
%                  velnorm_x = smooth(Tracks{itrack}(:,4),10);
                velnorm_z = smoothdata(Tracks{itrack}(:,3), 'movmean', 10);
                velnorm_x = smoothdata(Tracks{itrack}(:,4), 'movmean', 10);
            end
            
            if strcmp(mode,'2D_vel_z')
                % encode the direction of the velocity in positive/negative value
                velnorm = velnorm.*sign(mean(Tracks{itrack}(:,3)));
            end

            % remove out of grid bubbles (ie. the grid is too small)
            outP = zeros(numel(pos_z_round),1);
            outP = or(outP,pos_z_round<1);outP = or(outP,pos_z_round>sizeOut(1));
            outP = or(outP,pos_x_round<1);outP = or(outP,pos_x_round>sizeOut(2));

            ind=(sub2ind(sizeOut,pos_z_round(~outP),pos_x_round(~outP)));
            % Tracks are counted only once per pixel, the unique keeps only 1 point per pixel
            ind = unique(ind);
            velnorm = velnorm(~outP);
            velnorm_z = velnorm_z(~outP);
            velnorm_x = velnorm_x(~outP);
            intBubble = intBubble(~outP);
            
            for j=1:size(ind,1)
                % increment the intensity count of the pixel crossed by a track by +1
                MatOut(ind(j))=MatOut(ind(j))+1;
                % The sum of velocities will be average with Matout at the end of the code.
                MatOut_vel(ind(j)) = MatOut_vel(ind(j))+velnorm(j);
                MatOut_vel_z(ind(j)) = MatOut_vel_z(ind(j))+velnorm_z(j);
                MatOut_vel_x(ind(j)) = MatOut_vel_x(ind(j))+velnorm_x(j);
                MatOut_int(ind(j)) = MatOut_int(ind(j)) + intBubble(j);
            end
            
        end

    otherwise
        error('Wrong mode selected')
end

if exist('MatOut_vel')
    % average velocity
    MatOut_vel(MatOut>0) =  MatOut_vel(MatOut>0)./MatOut(MatOut>0);
    MatOut_vel_z(MatOut>0) =  MatOut_vel_z(MatOut>0)./MatOut(MatOut>0);
    MatOut_vel_x(MatOut>0) =  MatOut_vel_x(MatOut>0)./MatOut(MatOut>0);
    MatOut_int(MatOut>0) = MatOut_int(MatOut>0)./MatOut(MatOut>0);
    if nargout>1
        varargout{1} = MatOut_vel;
        varargout{2} = MatOut_vel_z;
        varargout{3} = MatOut_vel_x;
        varargout{4} = MatOut_int;
    else
        MatOut = MatOut_vel;
    end
end

end
