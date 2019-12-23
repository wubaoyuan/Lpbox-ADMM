function ha = shBreakLn(Lns, boxG, varargin)
% Show breaking line.
%
% Input
%   Lns      -  breaking line, 1 x (m - 1) (cell), 2 x 2
%   boxG     -  global bounding box, 2 x 2
%   varargin
%     show option
%     parAx  -  axis parameter, {[]} | ...
%
% Output
%   ha       -  handle
%     hLn    -  trajectory, 1 x k (cell), 2 x nFi
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 01-07-2012

% show option
psSh(varargin);

% function option
parAx = ps(varargin, 'parAx', []);

hold on;

% plot breaking line
ha.hLn = plot(Lns{1}(1, :), Lns{1}(2, :), '--k');

% axis
box = xBox(boxG, parAx);
setAx(box, parAx);
