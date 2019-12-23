function h = shLnIni(x, y, varargin)
% Show line.
%
% Input
%   x        -  endpoint position in x, 1 x 2
%   y        -  endpoint position in y, 1 x 2
%   varargin
%     show option
%     cl     -  line color, {'y'}
%     lnWid  -  line width, {1}
%
% Output
%   h        -  figure content handle
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 02-03-2010
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% show option
psSh(varargin);

% function option
cl = ps(varargin, 'cl', 'y');
lnWid = ps(varargin, 'lnWid', 1);

h = plot(x, y, '-', 'Color', cl, 'LineWidth', lnWid);
