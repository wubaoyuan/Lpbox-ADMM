function crbf_plot_grid(x, z, w, sigma_kernel)
% Plot RBF warpped grid.
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-03-2012

resolution = 4;
resolution_grid = 3;

[grid_pts, controls] = ctps_plot_grid_gen(x, resolution, resolution_grid);

% warp grid
grid_new = crbf_warp_pts(grid_pts, z, w, sigma_kernel);

% plot
ori_color = ones(1, 3) * 0.7;
new_color = 'b';

ctps_plot_gridbox(1, grid_pts, controls, ori_color, ':');
hold on;

ctps_plot_gridbox(1, grid_new, controls, new_color, '-');
hold on;
