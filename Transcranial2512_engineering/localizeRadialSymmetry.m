function [zc,xc] = localizeRadialSymmetry(I,fwhmz,fwhmx)
%% function [zc,xc] = localizeRadialSymmetry(I,fwhmz,fwhmx)
% Performs localization using radial symmetry properties
%
% Created by Baptiste Heiles on 05/09/18
% Inspired from Raghuveer Parthasarathy, The University of Oregon
%
% DATE 2020.07.22 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code published in 2020 article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.
% Open Platform for Ultrasound Localization Microscopy: performance assessment of localization algorithms
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018
%
% Calculates the center of a 2D intensity distribution.
% Method: The gradient of a function that has perfect radial symmetry will
% point towards the origin. Thus we take the local gradient and construct
% lines through any point with orientation parallel to the local gradient.
% The origin is the point that will minimize the distance between itself
% and all such lines.
% 具有完美径向对称性的函数的梯度将指向原点。
% 因此，我们取局部梯度，并通过任意一点构造方向与局部梯度平行的直线。
% 原点是能将自身与所有此类直线之间距离最小化的点。
% INPUTS:
%       - I : 2D intensity distribution
%        Size need not be an odd number of pixels along each dimension
%       - fwhmz, fwhmx : full width at half maximum in direction z and x (unused)
% OUTPUTS:
%       - [zc, xc] : the center of radial symmetry,
%            px, from center

%% Number of grid points
[Nz,Nx] = size(I);

%% Radial symmetry algorithm

% grid coordinates are -n:n, where Nz (or Nx) = 2*n+1
% grid midpoint coordinates are -n+0.5:n-0.5. Note that z increases "downward"
zm_onerow = (-(Nz-1)/2.0+0.5:(Nz-1)/2.0-0.5)';
zm = zm_onerow(:,ones(Nx-1, 1));
xm_onecol = (-(Nx-1)/2.0+0.5:(Nx-1)/2.0-0.5);
xm = xm_onecol(ones(Nz-1, 1),:);

% Calculate derivatives along 45-degree shifted coordinates (u and v) Please refer to Appendix 2 of the publication attached to this code for basis definition
dIdu = I(1:Nz-1,2:Nx)-I(2:Nz,1:Nx-1);% Gradient along the u vector
dIdv = I(1:Nz-1,1:Nx-1)-I(2:Nz,2:Nx);% Gradient along the v vector

% Smoothing the gradient of the I window
h = ones(3)/9;% 均值滤波器
fdu = conv2(dIdu, h, 'same');% Convolution of the gradient with a simple averaging filter
fdv = conv2(dIdv, h, 'same');% 二维卷积，把窗口内 9 个像素的值加起来除以 9
dImag2 = fdu.*fdu + fdv.*fdv; % Squared gradient magnitude 勾股定理

% Slope of the gradient . Please refer to appendix 2 of the publication attached to this code for basis/orientation
m = -(fdv + fdu) ./ (fdu-fdv);% 斜率

% Check if m is NaN (which can happen when fdu=fdv). In this case, replace with the un-smoothed gradient.
NNanm = sum(isnan(m(:)));
if NNanm > 0
    unsmoothm = (dIdv + dIdu) ./ (dIdu-dIdv);
    m(isnan(m))=unsmoothm(isnan(m));%大部分地方是高质量的平滑斜率，个别坏点用原始斜率顶替
end

% If it's still NaN, replace with zero and we'll deal with this later
NNanm = sum(isnan(m(:)));
if NNanm > 0
    m(isnan(m))=0;
end

% Check if m is inf (which can happen when fdu=fdv).
try
    m(isinf(m))=10*max(m(~isinf(m)));%找出当前所有正常斜率里的最大值，然后乘以 10
catch
    % Replace m with the unsmoothed gradient
    m = (dIdv + dIdu) ./ (dIdu-dIdv);
end

% Calculate the z intercept of the line of slope m that goes through each grid midpoint
b = zm - m.*xm; % y = kx+b

% Weight the intensity by square of gradient magnitude and inverse
% distance to gradient intensity centroid. This will increase the intensity of areas close to the initial guess
% 找这 16 条线交汇最密集的那个点。那个点就是微泡的亚像素中心。
sdI2 = sum(dImag2(:));
zcentroid = sum(sum(dImag2.*zm))/sdI2;% Initial guess of the centroid in z纵向重心
xcentroid = sum(sum(dImag2.*xm))/sdI2;% Initial guess of the centroid in x横向重心
w  = dImag2./sqrt((zm-zcentroid).*(zm-zcentroid)+(xm-xcentroid).*(xm-xcentroid));% 权重 距离远近/梯度强度

% least-squares minimization to determine the translated coordinate
% system origin (xc, yc) such that lines y = mx+b have
% the minimal total distance^2 to the origin:
% See function lsradialcenterfit (below)
[zc,xc] = lsradialcenterfit(m, b, w);

end

% We'll code the least square solution function separately as we could find the solution with another implementation
function [zc,xc] = lsradialcenterfit(m, b, w)
    % least squares solution to determine the radial symmetry center

    % inputs m, b, w are defined on a grid
    % w are the weights for each point
    wm2p1 = w./(m.*m+1); % 几何缩放
    sw  = sum(sum(wm2p1));
    smmw = sum(sum(m.*m.*wm2p1)); % 斜率平方的加权总和
    smw  = sum(sum(m.*wm2p1)); % 斜率的加权总和
    smbw = sum(sum(m.*b.*wm2p1));% 斜率与截距乘积的加权总和
    sbw  = sum(sum(b.*wm2p1));% 截距的加权总和
    det = smw*smw - smmw*sw;% 计算系数矩阵的行列式
    xc = (smbw*sw - smw*sbw)/det;    % relative to image center
    zc = (smbw*smw - smmw*sbw)/det; % relative to image center

end
