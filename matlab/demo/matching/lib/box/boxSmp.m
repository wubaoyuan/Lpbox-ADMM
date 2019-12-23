function X = boxSmp(box, n)
% Sampling point within the box.
% 
% Input
%   box     -  box, 2 x 2
%   n       -  #points
%
% Output
%   X       -  points, 2 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-13-2012

mi = box(:, 1);
ma = box(:, 2);
siz = ma - mi;
d = size(box, 1);

X = rand(d, n);
X = X .* repmat(siz, 1, n) + repmat(mi, 1, n);