function wsPath = cmumPath(src)
% Obtain file paths for the given CMU Motion source.
%
% Input
%   src     -  cmum src
%
% Output
%   wsPath
%     vdo   -  vdo path
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-04-2013

global footpath; % specified in addPath.m
foldpath = sprintf('%s/data/cmum', footpath);
wsPath.vdo = sprintf('%s/%s/images', foldpath, src.tag);
