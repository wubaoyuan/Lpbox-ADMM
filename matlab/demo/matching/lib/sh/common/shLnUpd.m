function shLnUpd(h, x, y)
% Show line.
%
% Input
%   h       -  figure content handle
%   x       -  endpoint position in x, 1 x 2
%   y       -  endpoint position in y, 1 x 2
%
% Output
%   h       -  figure content handle
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-03-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

set(h, 'XData', x, 'YData', y);
