function blocks = vec2block(a, n)
% Divide the vector into almost length-equal blocks.
%
% Input
%   a       -  the vector, 1 x na
%   n       -  #blocks
%
% Output
%   blocks  -  the blocks, 1 x n (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-21-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

na = length(a);
w = floor(na / n);

blocks = cell(1, n);
j = 0;
for i = 1 : n
    wi = select(i == 1, na - (w * (n - 1)), w);
    blocks{i} = a(j + 1 : j + wi);
    j = j + wi;
end
