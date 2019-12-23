function [angs, mags] = derX(X)
% Calculate the derative of curve in 2-D.
%
% Input
%   X       -  position, 2 x n
%
% Output
%   angs    -  direction of derivative, 1 x n
%   mags    -  magtitude of derivative, 1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 05-27-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

XD = gradient(X);
angs = atan2(XD(2, :), XD(1, :)); 
angs = angs / pi * 180;
mags = sqrt(sum(XD .^ 2));
