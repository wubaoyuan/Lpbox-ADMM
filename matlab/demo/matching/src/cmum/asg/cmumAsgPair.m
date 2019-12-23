function [nm, gaps, PFs, nIns] = cmumAsgPair(tag)
% Obtain all pairs of CMU Motion source.
%
% Input
%   tag     -  type of pair, 1 | 2 | 3 | 4
%
% Output
%   nm      -  'house' | 'hotel'
%   gaps    -  gap of baseline, 1 x mGap
%   PFs     -  frame index of each pair, 1 x mGap (cell), 2 x nP
%   nIns    -  #inlier, 1 x 2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-04-2013

if tag == 1
    gaps = 0 : 10 : 90;
    nIn = 30;
    nm = 'house';
    n = 111;

elseif tag == 2
    gaps = 0 : 10 : 90;
    nIn = 25;
    nm = 'house';
    n = 111;
    
elseif tag == 3
    gaps = 0 : 10 : 90;
    nIn = 30;
    nm = 'hotel';
    n = 101;

elseif tag == 4
    gaps = 0 : 30 : 90;
    nIn = 25;
    nm = 'hotel';
    n = 101;

else
    error('unknown tag: %d', tag);
end
mGap = length(gaps);
nIns = zeros(1, 2) + nIn;

% per gap
PFs = cell(1, mGap);
A = ones(n, n);
for iGap = 1 : mGap
    gap = gaps(iGap);

    Vis = triu(A, gap) - triu(A, gap + 1);
    [is, js] = find(Vis);
    PFs{iGap} = [is, js]';
end
