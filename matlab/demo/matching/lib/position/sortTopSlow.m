function [bs, idx] = sortTopSlow(as, m)
% Pick the top m elements.
%
% Input
%   as      -  all elements, 1 x n
%   m       -  #items
%
% Output
%   idx     -  integer, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-22-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 06-25-2012

[bs, idx] = zeross(1, m);
mi = min(as);

for i = 1 : m
    [b, p] = max(as);
    
    bs(i) = b;
    idx(i) = p;
    as(p) = mi;
end
