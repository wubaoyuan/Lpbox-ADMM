function [ids, idx] = uniId(id0s, m)
% Pick the top m id that are uniques.
%
% Input
%   id0s    -  all elements, 1 x n
%   m       -  #items
%
% Output
%   ids     -  selected elements, 1 x m
%   idx     -  #index of the elements
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-22-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-31-2012

% dimension
n = length(id0s);

[ids, idx] = zeross(1, m);
p = 1;

for i = 1 : m
    
    while true
    % check history
        isExist = 0;
        for j = 1 : i - 1
            if ids(j) == id0s(p)
                isExist = 1;
                break;
            end
        end
        
        if ~isExist
            break;
        end
        
        p = p + 1;
    end
    
    % store
    ids(i) = id0s(p);
    idx(i) = p;
    p = p + 1;
end
