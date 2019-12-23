function ctps_plot_gridbox(grid_method, grid, controls, marker_color, marker_type)
% Plot the grid box for TPS demonstration.
%
% Input
%   grid_method: 0 -- pts
%                1 -- pts + linking lines.
%   controls:    rows, cols, points_row, points_col
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 11-20-2012

mkSiz = 10;
    
rows = controls(1);
cols = controls(2);
points_row = controls(3);
points_col = controls(4);

if grid_method == 0
    plot(grid(:, 1), grid(:, 2), '.', 'color', marker_color, 'markersize', 5);
    
% used
elseif grid_method == 1
    for i = 1 : rows
        tmp = grid((i - 1) * points_row + 1 : i * points_row, :);
        ha = plot(tmp(:, 1), tmp(:, 2), 'color', marker_color, 'linestyle', marker_type);
        hold on;
        
        set(ha, 'markersize', mkSiz);
    end

    start_index = rows * points_row;
    for j = 1 : cols
        tmp = grid((j - 1) * points_col + start_index + 1 : j * points_col + start_index, :);
        ha = plot(tmp(:, 1), tmp(:, 2), 'color', marker_color, 'linestyle', marker_type);
        hold on;
        
        set(ha, 'markersize', mkSiz);
    end

else
    disp ('ERROR: cplot_gridbox -- wrong input parameters');
end
