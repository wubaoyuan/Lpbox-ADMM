function idx = idInGroup(ids, id0s)
% Find the index of a smaller group ("ids") in the larger group ("id0s").
%
% Input
%   ids     -  smaller group of ids, 1 x m
%   id0s    -  larger group of ids, 1 x m0
%
% Output
%   idx     -  #index of the elements, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 04-05-2013
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-05-2013

% dimension
m = length(ids);
m0 = length(id0s);

D = conDst(ids, id0s);
ind = find(D == 0);
[~, idx] = ind2sub(size(D), ind);