function [ids, idLefts, ms] = cluIdTop(id0s, G, kTop)
% Pick the top m id that are uniques.
%
% Input
%   id0s    -  all elements, 1 x n
%   G       -  label, k x n
%   kTop    -  #top cluster
%
% Output
%   ids     -  selected elements, 1 x m
%   idx     -  #index of the elements
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-22-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-02-2012

% dimension
n = length(id0s);

GTr = G(:, id0s);
ms = sum(GTr, 2);

% pick the top k
[~, idxClu] = sort(ms, 'descend');
GTr = GTr(idxClu(1 : kTop), :);
vis = sum(GTr, 1) > 0;

ids = id0s(vis);
idLefts = id0s(~vis);
