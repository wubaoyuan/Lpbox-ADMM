function sh3d(X, varargin)
% Show points in 3-D.
%
% Input
%   X        -  sample matrix, dim x n
%   varargin
%     G      -  class indicator matrix, {[]} | k x n
%     itIdx  -  index of interesting points, {[]}
%     mks    -  predefined mks, {[]}
%     mkSiz  -  size of mks, {5}
%     leg    -  legends, {[]}
%     anotI  -  index of point, {[]}
%     anotS  -  string of annotation, {[]}
%     face   -  flag of showing face color, 'y' | {'n'}
%
% Output
%   h        -  figure content handle
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 12-30-2008
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% show option
psSh(varargin);

% function option
G = ps(varargin, 'G', []);
mk = ps(varargin, 'mk', 'o');
cl = ps(varargin, 'cl', 'r');
mkSiz = ps(varargin, 'mkSiz', 5);
leg = ps(varargin, 'leg', []);
anotI = ps(varargin, 'anotI', []);
anotS = ps(varargin, 'anotS', []);
isFace = psY(varargin, 'face', 'y');
itIdx = ps(varargin, 'itIdx', []);

% dimension
n = size(X, 2);
X = pca(X, struct('dim', 3, 'debg', 'n'));

% label
if isempty(G)
    G = ones(1, n); 
end
k = size(G, 1);

% marker & color
if length(mk) == 1
    mk = cellRep(mk, k);
end
if length(cl) == 1
    cl = cellRep(cl, k);
end
if length(mkSiz) == 1
    mkSiz = repmat(mkSiz, 1, k);
end

% main plot
hold on;
pois = cell(1, k);
for c = 1 : k
    vis = G(c, :) == 1;

    pois{c} = plot3(X(1, vis), X(2, vis), X(3, vis), mk{c}, 'Color', cl{c}, 'MarkerSize', mkSiz(c), 'MarkerEdgeColor', cl{c});
    
    if isFace
        set(pois{c}, 'MarkerFaceColor', cl{c});
    end
end

% boundary
% setBound(X, 'mar', [.1 .1 .1]);

% angle
% view([15 20]);
