function m = rpmCorOpt(vx, y, T, algX, m_outliers_row, m_outliers_col, it_total, sigma)
% Update correspondence for RPM.
%
% Input
%   vx        -  data points, xmax x 2
%   y         -  data points, ymax x 2
%   T         -  temperature
%   algX      -  method name, 'icp' | 'mixture' | 'mix-rpm' | 'rpm' | 'rpm-old'
%   m_outliers_row
%   m_outliers_col
%   it_total
%   sigma
%
% Output
%   m         -  soft correspondence, xmax x ymax
%
% History
%   create    -  Anand Rangarajan (anand@noodle.med.yale.edu), 04-27-2000
%   modify    -  Feng Zhou (zhfe99@gmail.com), 04-23-2012

% dimension
[xmax, dim] = size(vx);
[ymax, dim] = size(y);

% ICP
if strcmp(algX, 'icp')
    m = rpmCorOptIcp(vx, y, sigma);
    m = m + randn(xmax, ymax) * (1 / xmax) * 0.001;

% one way mixture
elseif strcmp(algX, 'mixture')

    % given v = tranformed(x), update m
    y_tmp = zeros(xmax, ymax);
    for it_dim = 1 : dim
        y_tmp = y_tmp + (vx(:, it_dim) * ones(1, ymax) - ones(xmax, 1) * y(:, it_dim)') .^ 2;
    end

    m_tmp = 1 / sqrt(T) .* exp(-y_tmp / T);
    m_tmp = m_tmp + randn(xmax, ymax) * (1 / xmax) * 0.001;

    m = m_tmp;

    % normalize accross the outliers as well:
    sy = sum(m) + m_outliers_row;
    m = m ./ (ones(xmax, 1) * sy);

% mixture - RPM
elseif strcmp(algX, 'mix-rpm')

    % Given v = tranformed(x), update m:
    y_tmp = zeros(xmax, ymax);
    for it_dim = 1 : dim
        y_tmp = y_tmp + (vx(:, it_dim) * ones(1, ymax) - ones(xmax, 1) * y(:, it_dim)') .^ 2;
    end

    m_tmp = 1 / sqrt(T) .* exp(-y_tmp / T);
    m_tmp = m_tmp + randn(xmax, ymax) * (1 / xmax) * 0.001;

    m = m_tmp;
    [m, junk1, junk2] = rpmCorNor(m_tmp, m_outliers_col, m_outliers_row);

% RPM, double normalization
elseif strcmp(algX, 'rpm')
    
    % Given v = tranformed(x), update m:
    y_tmp = zeros(xmax, ymax);
    for it_dim = 1 : dim
        y_tmp = y_tmp + (vx(:, it_dim) * ones(1, ymax) - ones(xmax, 1) * y(:, it_dim)') .^ 2;
    end

    m_tmp = exp(-y_tmp / T);
    m_tmp = m_tmp + randn(xmax, ymax) * (1 / xmax) * 0.001;
    mi = min(m_tmp(:));
    if mi < -eps
        m_tmp = m_tmp - mi;
    end
    
    % double normalization, but keep outlier entries constant.
    moutlier = 1 / xmax * 0.1;
    m_outliers_row = ones(1, ymax) * moutlier;
    m_outliers_col = ones(xmax, 1) * moutlier;

    [m, junk1, junk2] = rpmCorNor(m_tmp, m_outliers_col, m_outliers_row);

else
    error('unknown correspondence method: %s', algX);
end
