function [L, Pt, rs, idx] = imgSeg(R0, k, varargin)
% Image segmentation and find the peak position in each segment.
%
% Input
%   R       -  response map, h x w
%   k       -  #clusters
%
% Output
%   L       -  label map, h x w
%   Pt      -  peak (center) of each cluster, 2 x k
%   rs      -  response value at each peak, 1 x k
%   idx     -  position of each peak, 1 x k
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-09-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-25-2012

% dimension
[h, w] = size(R0);

% re-scale
mi = min(R0(:));
ma = max(R0(:));
R = (R0 - mi) / (ma - mi) * max(h, w) * 3;

% position
xs = 1 : w;
ys = 1 : h;
[X, Y] = meshgrid(xs, ys);

% feature
Z = [X(:), Y(:), R(:)]';

% clustering
[~, labs] = kmeanFast(Z', k);

% label map
L = reshape(labs, [h, w]);

% per cluster
Pt = zeros(2, k);
[rs, idx] = zeross(1, k);
for c = 1 : k
    idxc = find(labs == c);

    % peak
    [~, p] = max(R(idxc));
    pc = idxc(p);
    
    % store
    Pt(:, c) = [X(pc); Y(pc)];
    rs(c) = R0(pc);
    idx(c) = pc;
end
