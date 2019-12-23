function [grid_pts, controls] = ctps_plot_grid_gen(x, res, res_grid)
% Generate grid points for displaying TPS deformation.
%
% Input
%   x       -  point set, n x 2
%   res     -  resolution, 4 | 
%   res_grid  -  grid, 3 | 
%
% Output
%   grid_pts  -  grid points, m x 2
%   controls  -  [rows; cols; points_row; points_col], 4 x 1
% 
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-21-2012

% set grid range
xrange = [min(x(:, 1)), max(x(:, 1))];
yrange = [min(x(:, 2)), max(x(:, 2))];

% expand a little bit
expand_ratio = 5;
xrange(1) = xrange(1) - (xrange(2) - xrange(1)) / expand_ratio;
xrange(2) = xrange(2) + (xrange(2) - xrange(1)) / expand_ratio;
yrange(1) = yrange(1) - (yrange(2) - yrange(1)) / expand_ratio;
yrange(2) = yrange(2) + (yrange(2) - yrange(1)) / expand_ratio;

% grid size
a = min(xrange(2) - xrange(1), yrange(2) - yrange(1)) / res;
grid_step = a / res_grid;

cols = ceil((xrange(2) - xrange(1)) / a + 1);
rows = ceil((yrange(2) - yrange(1)) / a + 1);

xrange(2) = xrange(1) + (cols - 1) * res_grid * grid_step;
yrange(2) = yrange(1) + (rows - 1) * res_grid * grid_step;

grid_pts = [];

% points_row, points_cols -- points along each row and col.
points_row = (cols - 1) * res_grid + 1; % two ending points.
for i = 1 : rows
    tmp_row = [[xrange(1) : grid_step : xrange(2)]', ...
               ones(points_row, 1) * (i - 1) * a + yrange(1)];
    grid_pts = [grid_pts; tmp_row];
end

points_col = (rows - 1) * res_grid + 1; % two ending points.
for j = 1 : cols
    tmp_col = [ones(points_col, 1) * (j - 1) * a + xrange(1), ...
               [yrange(1) : grid_step : yrange(2)]'];
    grid_pts = [grid_pts; tmp_col];
end

% store
controls = zeros(4, 1);
controls(1) = rows;
controls(2) = cols;
controls(3) = points_row;
controls(4) = points_col;
