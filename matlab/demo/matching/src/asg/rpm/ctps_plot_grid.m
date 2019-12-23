function ctps_plot_grid(x, y, c, d)
% Display TPS deformed grid.
%
% Input
%   x       -  TPS basis points.
%   y       -  points to be warped.
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-20-2012

resolution = 4;
resolution_grid = 3;

% generate grid points
[grid_pts, controls] = ctps_plot_grid_gen(x, resolution, resolution_grid);

% warp the grid
grid_new = ctps_warp_pts(grid_pts, x, c, d);

% plot
ori_color = ones(1, 3) * 0.7;
new_color = 'b';

% plot the original grid
ctps_plot_gridbox(1, grid_pts, controls, ori_color, ':'); 
hold on;

% plot the new grid
ctps_plot_gridbox(1, grid_new, controls, new_color, '-');
hold on;
