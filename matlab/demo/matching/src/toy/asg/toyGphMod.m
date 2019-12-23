function [mes, Vars] = toyGphMod(shp, k)
% Generate toy graph models.
%
% Input
%   shp     -  graph shape, 1 | ...
%              1: random graph
%   k       -  #nodes
%
% Output
%   mes     -  mean of Gaussian, 1 x k (cell), 2 x 1
%   Vars    -  variance of Gaussian, 1 x k (cell), 2 x 2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-16-2012

% random graph
if shp == 1
    % Gauss for each node
    [mes, Vars] = modGausses(k, 'dis', 'unif', 'maMiD', 'n', ...
                             'rBd', [.007, .02]);

else
    error('unknown shape: %d', shp);
end
