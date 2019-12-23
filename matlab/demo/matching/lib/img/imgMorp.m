function F = imgMorp(F0, par)
% Image morphological operation.
%
% Input
%   F0      -  original image
%   par     -  parameter
%     morp  -  radius, [] | 1 | 2 | ...
% 
% Output
%   F       -  new image
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-01-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function parameter
morp = ps(par, 'morp', []);

if isempty(morp)
    F = F0;
else
    se = strel('disk', morp);
    F1 = imerode(F0, se);
    F = imdilate(F1, se);
end
