function m = rpmCorOptIcp(vx, y, sigma)
% Calculate the correspondence of ICP.
%
% Input
%   vx      -  1st point set
%   y       -  2nd point set
%   sigma   -  sigma
%
% Output
%   m       -  soft correspondence, xmax x ymax
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-01-2012

% dimension
[xmax, dim] = size(vx);
ymax = size(y, 1);

% find nearest neighbour for each pt in x:
[M1, dist_x] = findneighbours(vx, y);
[M2, dist_y] = findneighbours(y, vx);

% set the distance threshold    
if sigma == 0
    distTh = 1e10;
else
    dist = [dist_x; dist_y];
    mean_x = mean(dist);
    sx = std(dist);
    distTh = mean_x + sigma * sx;
end

% correspondence
[m1, m2] = zeross(xmax, ymax); 
for i = 1 : xmax
    if dist_x(i) <= distTh
        m1(i, M1(i)) = 1;
    end
end
for j = 1 : ymax
    if dist_y(j) <= distTh
        m2(M2(j), j) = 1; 
    end
end
m = (m1 + m2) / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M, dists] = findneighbours(x, t)
% Find nearest neighbour in template t for each point in x.
%
% Input
%   x      -  data points, m x dim
%   t      -  data points, n x dim
%
% Output   
%   M      -  index, m x 1
%             M(1) = 10, y(10) is nearest from x(1);
%   dists  -  distance, m x 1

% dimension
[m, dim] = size(x);
[n, dim] = size(t);
M = zeros(m, 1);

% |x - t| matrices
xttmp = zeros(n, m);
for i = 1 : dim
    xtmp = ones(n, 1) * x(:, i)';
    ttmp = t(:, i)  * ones(1, m); 
    xttmp = xttmp + (xtmp - ttmp) .* (xtmp - ttmp);
end

% M + min_dist list
[min_dist, min_index] = min(xttmp);
dists = (sqrt(min_dist))';
M = min_index';
