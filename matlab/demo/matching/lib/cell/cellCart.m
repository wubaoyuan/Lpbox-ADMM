function varargout = cellCart(pos, varargin)
% Return the Cartesian prodcut of sets.
%
% Example
%   input      -  a0s = {1, 2}; b0s = {'aa', 'bb'}; c0s = {'name', 'sex'};
%   call       -  [as, bs, cs] = cellCart([1 1 2], a0s, b0s, c0s)
%   output     -  as = {1, 2, 1, 2};
%              -  bs = {'aa', 'bb', 'aa', 'bb'};
%              -  cs = {'name', 'name', 'sex', 'sex'};
%
% Input
%   pos        -  position, 1 x m
%   varargin   -  set of cell array, 1 x m (cell)
%
% Output
%   varargout  -  set of cell array, 1 x m (cell), 1 x n (cell)
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 02-16-2009
%   modify     -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% position
m = length(varargin);
if isempty(pos)
    pos = 1 : m;
end
k = max(pos);
[vis, ns] = zeross(1, k);

n0s = cellDim(varargin, 2);
n = 1;
for j = 1 : m
    pj = pos(j);
    nj = n0s(j);

    if vis(pj) == 0
        vis(pj) = 1;
        ns(pj) = nj;
        n = n * nj;
    else
        if ns(pj) ~= nj
            error('inconsistent dimension');
        end
    end
end

for j = 1 : m
    varargout{j} = cell(1, n);
end

subs = cell(1, m);
for i = 1 : n
    [subs{:}] = ind2sub(ns, i);

    for j = 1 : m
        pj = pos(j);
        varargout{j}{i} = varargin{j}{subs{pj}};
    end
end
