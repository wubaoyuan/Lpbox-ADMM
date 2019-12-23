function idx = finds(a, b)
% Find the position of b's element in a.
%
% Sample
%   input   -  a = [3, 5, 1, 4, 2]; b = [1 2]
%   call    -  idx = finds(a, b)
%   output  -  b = [3 5];
%
% Input
%   a       -  vector a, 1 x n
%   b       -  vector b, 1 x m
%
% Output
%   idx     -  position, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-03-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

m = length(b);

idx = zeros(1, m);
for i = 1 : m
    tmp = find(a == b(i));
    if ~isempty(tmp)
        idx(i) = tmp;
    end
end
