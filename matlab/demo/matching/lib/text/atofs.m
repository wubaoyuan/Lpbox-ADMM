function fs = atofs(as)
% Convert multiple strings to float array.
%
% Input
%   as      -  string set, 1 x m (cell)
%
% Output
%   fs      -  double array, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-25-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m = length(as);

fs = zeros(1, m);
for i = 1 : m
    fs(i) = atof(as{i});
end
