function [head, colMa, Col] = imgCol(FT, F)
% Normalized Cross-correlation.
%
% Input
%   FT      -  template image, hT x wT
%   F       -  source image, h x w
%
% Output
%   head    -  upper-left corner of the part in source image, 2 x 1
%   colMa   -  max value of cross-correlation
%   Col     -  value of cross-correlation, h x w
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
sizT = size(FT); 
siz = size(F);

% normalized cross-correlation
Col = normxcorr2(FT, F);

marL = floor((sizT + 1) / 2);
marR = sizT - marL;

Col = Col(marL(1) : marL(1) + siz(1) - 1, marL(2) : marL(2) + siz(2) - 1);
[colMa, idx] = max(Col(:));

% upper-left corner
[i, j] = ind2sub(size(Col), idx(1));
head = [i, j] - marR + 1;
head = head';
