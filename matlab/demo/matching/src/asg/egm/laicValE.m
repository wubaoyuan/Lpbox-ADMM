function val = laicValE(E)
% Make sure every point has three neighbors.
%
% Input
%   E       -  edge connection, n x n
%
% Output
%   val     -  flag of whether point is valid, 1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-25-2012

% dimension
n = size(E, 1);
val = ones(1, n);

ks = zeros(1, n);
for nIt = 1 : n
    for i = 1 : n
        if val(i) == 0
            continue;
        end
        ks(i) = length(find(E(i, :)));
    end

    [kMi, p] = min(ks);
    
    % remove this point
    if kMi < 3
        val(p) = 0;
        ks(p) = n;
    else
        break;
    end
end

val = val == 1;