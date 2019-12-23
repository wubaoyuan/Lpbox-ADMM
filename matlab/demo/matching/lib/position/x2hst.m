function hst = x2hst(x, s)
% Obtain histogram.
%
% Input
%   x       -  a set of instances, 1 x n
%   s       -  starting position of each segment, 1 x (m + 1)
%
% Output
%   hst     -  histogram, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-30-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-06-2013

% dimension
n = length(x);
m = length(s) - 1;

hst = zeros(1, m);
for i = 1 : m
    if i < m
        hst(i) = length(find(x >= s(i) & x < s(i + 1)));
    else
        hst(i) = length(find(x >= s(i) & x <= s(i + 1)));
    end
end
