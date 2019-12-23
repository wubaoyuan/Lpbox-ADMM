function img = imgBlur(img0, par)
% Image blurring.
%
% Input
%   img0    -  initial image
%   par     -  parameter
%     blur  -  radius of gaussian filter, {[]} | 1 | 2 | ...
% 
% Output

%   img     -  transformed image
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-24-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function parameter
blur = ps(par, 'blur', []);

if isempty(blur)
    img = img0;

else
    H = fspecial('gaussian', blur * 4, blur);
    img = imfilter(img0, H, 'replicate');
end
