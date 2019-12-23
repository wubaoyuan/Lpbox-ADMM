function Rot = a2rot(a)
% Obtain the 2-D rotation matrix from the given angle.
%
% Input
%   a       -  angle
%
% Output
%   Rot     -  rotation matrix, 2 x 2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-19-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

Rot = [cos(a), sin(a); ...
      -sin(a), cos(a)];

