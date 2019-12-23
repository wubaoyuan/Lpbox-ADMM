function p = inSt(iF, R)
% Obtain the position in the specified range of frame.
%
% Input
%   iF      -  frame position
%   R       -  range matrix, 2 x m
%
% Output
%   p       -  range position
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

m = size(R, 2);
for i = 1 : m
    if R(1, i) <= iF && iF <= R(2, i)
        p = i;
        return;
    end
end
p = 0;
