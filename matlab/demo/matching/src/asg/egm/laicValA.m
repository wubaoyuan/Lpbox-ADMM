function [A, mR] = laicValA(A0)
% Make sure every point has three neighbors.
%
% Input
%   A0      -  original edge connection, n x n
%
% Output
%   A       -  new edge connection, n x n
%   mR      -  #rows that have been updated
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-27-2012

% dimension
n = size(A0, 1);

A = A0;
mR = 0;
for i = 1 : n
    vis = A0(i, :) == 1;
    m = length(find(vis));
    
    if m >= 3
        continue;
    end
    
    while true
        idx = randperm(n);
        idx = idx(1 : 3 - m);
        
        if all(~vis(idx))
            break;
        end
    end
    
    A(i, idx) = 1;
    A(idx, i) = 1;
    mR = mR + 1;
end

if mR > 0
    pr('%d rows has been updated', mR);
end
