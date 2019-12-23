function X = boxVex(box)
% Obtain all vertices of a given box.
%
% Input
%   box     -  box, dim (= 3) x 2
%
% Output
%   X       -  box vertices, dim (= 3) x 8
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 11-05-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

X = zeros(3, 8);
ii = 0;
for i = 1 : 2
    for j = 1 : 2
        for k = 1 : 2
            ii = ii + 1;
            X(1, ii) = box(1, i);
            X(2, ii) = box(2, j);
            X(3, ii) = box(3, k);
        end
    end
end
