function [bas, ba0s] = laicHull(Pt, C)
% Core implementation of LAIC algorithm.
%
% Input
%   Pt      -  point set, 2 x n
%   C       -  feature matching cost matrix, k x n
%
% Output
%   bas     -  basis, 1 x k (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-28-2012

% dimension
[k, n] = size(C);

% initial basis
[bas, ba0s] = cellss(k, 1);
for i = 1 : k
%     bas{i} = 1 : n;
% 
%     % convex hull
%     costSurf = [Pt(bas{i}, :), C(i, bas{i})'];
%     costSurfHull = convhulln(costSurf, {'Qt' 'PD2'});
% 
%     % support points of the convex hull
%     costSurfHull = sort(costSurfHull(:));
%     ind = costSurfHull ~= [costSurfHull(2 : end); 0];
%     bas{i} = costSurfHull(ind);

    % visualize by feng
%     figure(10); clf;
%     trisurf(costSurfHull, costSurf(:, 1), costSurf(:, 2), costSurf(:, 3));
%     hold on;
%     plot3(costSurf(:, 1), costSurf(:, 2), costSurf(:, 3), 'ro', 'MarkerSize', 10);
%     plot3(costSurf(bas{1}, 1), costSurf(bas{1}, 2), costSurf(bas{1}, 3), 'b+', 'MarkerSize', 10);
    bas{i} = find(C(i, :) < 1e4);
    ba0s{i} = 1 : n;
end
