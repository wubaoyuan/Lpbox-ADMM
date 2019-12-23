function [m, m_outlier_col, m_outlier_row] = rpmCorNor(m, m_outlier_col, m_outlier_row)
% Double normalization of m (with outlier col/row).
%
% Input
%   m
%   m_outlier_col
%   m_outlier_row
%
% Output
%   m
%   m_outlier_col
%   m_outlier_row
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-23-2012

% parameters
norm_threshold = 0.05;
norm_maxit = 10;

% dimension
[xmax, ymax] = size(m);

norm_it = 0;
while norm_it <= norm_maxit
    % row normalization
    sx = sum(m')' + m_outlier_col;
    m  = m ./ (sx * ones(1, ymax));
    m_outlier_col = m_outlier_col ./ sx;
  
    % column normalization
    sy = sum(m) + m_outlier_row;
    m  = m ./ (ones(xmax, 1) * sy);
    m_outlier_row = m_outlier_row ./ sy;
  
    % time to quit?
    err = ((sx - 1)' * (sx - 1) + (sy - 1) * (sy - 1)') / (length(sx) + length(sy));
    if err < (norm_threshold .* norm_threshold)
        break;
    end

    % run out of time
    norm_it = norm_it + 1;
end

% mi = min(m(:));
% if mi < -1e-3
%     fprintf('error %.3f\n', mi);
% end
