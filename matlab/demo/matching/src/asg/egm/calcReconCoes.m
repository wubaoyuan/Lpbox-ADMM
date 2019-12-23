function [H, weis, Err] = calcReconCoes(Pts, As, Use, eta)
% Calculates the reconstruction weights for each point given a point-set.
%
% Input
%   Pts     -  point set, 1 x m (cell), d x k
%   As      -  graph adjacency matrix, 1 x m (cell), k x k (binary)
%   Use     -  flag of used graphs, k x m
%   eta     -  eta for reguralization, {1000}
%
% Output    
%   H       -  reconstruction matrix, H = I - W, k x k
%   weis    -  weights, 1 x k
%   Err     -  error, k x m
%
% History   
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-04-2012

% dimension
m = length(Pts);
[d, k] = size(Pts{1});
prIn('calcReconCoes', 'm %d, k %d, d %d', m, k, d);

% neighbors
Vis = zeros(k, k);
for i = 1 : m
    Vis = Vis | As{i} ~= 0;
end
ns = sum(Vis, 2);
idxs = cell(1, k);
for c = 1 : k
    idxs{c} = find(Vis(c, :));
end

% normalize the scale
Pt0s = Pts;
Pts = falNormSca(Pt0s, []);

% calculate the convex combination for each point
H = zeros(k, k);
for c = 1 : k
    vis = Use(c, :) == 1;
    H(c, :) = calcReconH(Pts(vis), idxs{c}, c, eta);
end

% error
Err = zeros(m, k);
for i = 1 : m
    Pt = Pts{i};
    PtD = Pt * H';

    Err(i, :) = sum(abs(PtD));
end
Err = Err';

% error -> weight (not very useful)
% errs = sum(Err, 2);
% errSum = sum(errs);
% weis = errSum - errs;
% weis = weis / sum(weis);
% weis = weis * 29;
weis = ones(1, 29);

prOut;
