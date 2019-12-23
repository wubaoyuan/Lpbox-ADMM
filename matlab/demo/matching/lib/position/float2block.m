function X = float2block(X0, k1, k2)
% Divide the [0, 1] into k segments. Round the elements of the given continuous matrix into nearest segments.
%
% Example
%    k = 3: 0.0 ~ 0.3 => 1
%           0.3 ~ 0.6 => 2
%           0.6 ~ 1.0 => 3
%
% Input
%   X0  -  continuous matrix with elements x in [0, 1]
%   k1  -  number of segments
%   k2  -  if indicated, mapping to [k1 k2]
%
% Ouput
%   X   -  discrete matrix
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-11-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if nargin == 2
    k2 = k1; k1 = 1;
end

if k2 < k1
    tmp = k1; k1 = k2; k2 = tmp;
end

k = k2 - k1 + 1;
X = floor(X0 * k + eps) + k1;
