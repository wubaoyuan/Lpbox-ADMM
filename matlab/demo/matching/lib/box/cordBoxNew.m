function cords = cordBoxNew(cord0s, box0s, boxs, Sca)
% Obtain new coordinate.
% 
% Input
%   cord0s  -  original coordinate, 1 x m (cell), 2 x kJ x nFi
%   box0s   -  original box set, 1 x m (cell), 2 x 2
%   boxs    -  new box set, 1 x m (cell), 2 x 2
%   Sca     -  scale (= siz / siz0), 2 x m
%
% Output
%   cords   -  new coordinate, 1 x m (cell), 2 x kJ x nFi
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-22-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m = length(cord0s);

% per instance
cords = cell(1, m);
for i = 1 : m
    cords{i} = xBoxNew(cord0s{i}, box0s{i}, boxs{i}, Sca(:, i));
end
