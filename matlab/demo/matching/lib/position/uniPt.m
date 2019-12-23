function [Pt, isReds] = uniPt(Pt0, varargin)
% Pick points that are unique.
%
% Input
%   Pt0     -  original points, 2 x n0
%   varargin
%     th    -  threshold, {1e-7}
%
% Output
%   Pt      -  new points, 2 x n
%   isReds  -  flag of redundance, 1 x n0
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-22-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-21-2012

% function option
th = ps(varargin, 'th', 1e-7);

% dimension
n0 = size(Pt0, 2);

% pair-wise distance
D = conDst(Pt0, Pt0);

% skip the diagonal
idx = sub2ind([n0 n0], 1 : n0, 1 : n0);
D(idx) = inf;

isReds = zeros(1, n0);
for i = 1 : n0
    if isReds(i)
        continue;
    end

    js = find(D(i, :) < th);
    
    if isempty(js)
        continue;
    end
    
    D(js, :) = inf;
    D(:, js) = inf;
    isReds(js) = 1;
end

Pt = Pt0(:, isReds == 0);