function shHst3d(Hist, varargin)
% Plot histogram in 3-D space.
%
% Input
%   Hist    -  histogram matrix, n x m
%   varargin
%     show option
%     zLim  -  threshold, {[]}
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-09-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
psShow(varargin);
zLim = ps(varargin, 'zLim', []);

[m, n] = size(Hist);
bar3(Hist');

% boundary
xlim([1, m] + [-.5, .5]);
ylim([1, n] + [-.5, .5]);
if isempty(zLim)
    mi = 0; ma = max(Hist(:));
else
    mi = zLim(1); ma = zLim(2);
end
zlim([mi, ma]);

view([30 30]);

set(gca, 'GridLineStyle', 'none');
set(gca, 'color', [1 1 1] * .5);
