function [idss, weis, ms] = cluId(id0s, G)
% Pick the top m id that are uniques.
%
% Input
%   id0s    -  all elements, 1 x n
%   G       -  label, k x n
%
% Output
%   idss    -  selected elements, 1 x k (cell), 1 x ni
%   weis    -  weights, 1 x k
%   ms      -  #elements, 1 x k
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-22-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-02-2012

% dimension
[k, n] = size(G);

GTr = G(:, id0s);
ms = sum(GTr, 2);

% sort by size
[ms, idxClu] = sort(ms, 'descend');
GTr = GTr(idxClu, :);

idss = cell(1, k);
for c = 1 : k
    idx = find(GTr(c, :));
    
    idss{c} = id0s(idx);
end

% weight
weis = ms / sum(ms);
