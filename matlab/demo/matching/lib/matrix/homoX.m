function Xs = homoX(Xs, flag)
% Making feature be homogeneous coordinates or removing the homogeneous 1s.
%
% Input
%   Xs      -  original sample matrix, 1 x m (cell), dim x n
%   flag    -  flag, 1 | 0
%              1: add one row of 1s at the bottom of each Xs{i}
%              0: remove the last row of each Xs{i}
%
% Output
%   Xs      -  new sample matrix, 1 x m (cell), (dim + 1) x n | (dim - 1) x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 11-19-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-20-2013

% default parameter
if ~exist('flag', 'var')
    flag = 1;
end

% dimension
m = length(Xs);

for i = 1 : m
    if flag == 1
        Xs{i} = [Xs{i}; ones(1, size(Xs{i}, 2))];
    else
        Xs{i}(end, :) = [];
    end
end
