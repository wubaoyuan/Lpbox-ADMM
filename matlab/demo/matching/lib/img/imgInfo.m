function [siz, nC, isU] = imgInfo(F)
% Obtain image information.
%
% Input
%   F       -  image, h x w x 3 (uint8) | h x w (uint8) | h x w (double)
% 
% Output
%   siz     -  size, [height width]
%   nC      -  #channels, 1 | 3
%   isU     -  flag of using uint8 format
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-02-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
if ndims(F) == 2
    [h, w] = size(F);
    nC = 1;
else
    [h, w, nC] = size(F);
end
siz = [h, w];

% data type
isU = isa(F, 'uint8');
