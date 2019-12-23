function F = imgNew(siz, nC, isU)
% Create an empty image.
% 
% Input
%   siz     -  size
%   nC      -  #channels, 1 | 3
%   isU     -  flag of using uint8 format
%
% Output
%   F       -  image, h x w x 3 (uint8) | h x w (uint8) | h x w (double)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-02-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if isU
    if nC == 1
        F = zeros(siz(1), siz(2), 'uint8');
    else
        F = zeros(siz(1), siz(2), nC, 'uint8');
    end
else
    if nC == 1
        F = zeros(siz(1), siz(2));
    else
        error('invalid image format');
    end
end
