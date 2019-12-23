function hist = x2hist(x, ma)
% Count the appearance number of the value stored in x.
% 
%
% Input
%   x       -  input vector, 1 x n
%
% Output
%   hist    -  count vector, 1 x max(x)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-15-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

n = length(x);

if ~exist('ma', 'var')
    ma = max(x);
end

hist = zeros(1, ma);
for i = 1 : n
    if x(i) == 0, continue; end
    hist(x(i)) = hist(x(i)) + 1;
end
