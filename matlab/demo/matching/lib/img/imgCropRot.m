function F = imgCropRot(F0, cen, rad, ang)
% Crop image with rotation.
%
% Input
%   F0      -  original image
%   cen     -  center point, 2 x 1
%   rad     -  radius
%   ang     -  rotation angle
%
% Output
%   F       -  cropped image
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-30-2012

% image size
[sizF0, nC, isU] = imgInfo(F0);

% original box pixel
x0s = -rad : rad;
y0s = -rad : rad;
[X0, Y0] = meshgrid(x0s, y0s);
Z0 = [X0(:), Y0(:)]';

% dimension
siz = size(X0);
n = siz(1) * siz(2);

% rotate the box
R = [cos(ang), sin(ang); -sin(ang), cos(ang)];
Z = R * Z0 + repmat(cen, 1, n);
X = reshape(Z(1, :), siz);
Y = reshape(Z(2, :), siz);

% interpolation
F = imgNew(siz, nC, isU);
for c = 1 : nC
    Fc = interp2(double(F0(:, :, c)), X, Y);
    
    if isU
        F(:, :, c) = uint8(Fc);
    else
        F(:, :, c) = Fc;
    end
end
