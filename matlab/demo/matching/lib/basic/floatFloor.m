function f = floatFloor(f0, precision)
% Floor for float numbers with the specified precision.
%
% Input
%   f0         -  input float number
%   precision  -  number of bit, e.g. 1 for 0.1, 2 for 0.01
%
% Output
%   f          -  output float number
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify     -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

tmp = 10 ^ precision;
f = floor((f0 + eps) * tmp) / tmp;
