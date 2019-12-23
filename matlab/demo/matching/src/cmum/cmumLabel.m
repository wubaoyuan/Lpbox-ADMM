function Ps = cmumLabel(tag)
% Load ground-truth label file for CMU Motion sequence.
%
% Input
%   tag     -  sequence name, 'house' | 'hotel'
%
% Output
%   Ps      -  label points, 1 x n (cell), 2 x 30
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-29-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-03-2013

% specified in addPath.m
global footpath;
foldpath = sprintf('%s/data/cmum', footpath);

% #images
if strcmp(tag, 'house')
    n = 111;

elseif strcmp(tag, 'hotel')    
    n = 101;

else
    error('unknown tag %s', tag);
end

% read all label
Ps = cell(1, n);
for i = 1 : n
    
    fpath = sprintf('%s/%s/label/%s%d', foldpath, tag, tag, i);

    P = load(fpath);
    
    % store
    Ps{i} = P;
end



