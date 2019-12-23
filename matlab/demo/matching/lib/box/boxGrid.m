function [X, Y, Z] = boxGrid(box, siz)
% Obtain 3-d grid for 3-d box.
%
% Input
%   box     -  box, dim (= 3) x 2
%   siz     -  number of bins, dim x 1
%
% Output
%   X       -  grid position in x
%   Y       -  grid position in y
%   Z       -  grid position in z
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-22-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

h = siz(1);
w = siz(2);

xs = linspace(box(1, 1), box(1, 2), w);
ys = linspace(box(2, 1), box(2, 2), w);
zs = linspace(box(3, 2), box(3, 1), h);
[X, ~] = meshgrid(xs, zs);
[Y, Z] = meshgrid(ys, zs);
