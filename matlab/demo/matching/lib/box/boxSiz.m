function Siz = boxSiz(boxs)
% Obtain the size of boxes.
% 
% Remark
%   Each box is a dim x 2 matrix, where 1st and 2nd column denote
%   the minimum and maximum value of each dimension respectively.
%
% Input
%   boxs    -  set of bounding boxs, 1 x m (cell), dim x 2
%
% Output
%   Siz     -  size of bounding boxs, dim x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-20-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m = length(boxs);
dim = size(boxs{1}, 1);

Siz = zeros(dim, m);
for i = 1 : m
    Siz(:, i) = boxs{i}(:, 2) - boxs{i}(:, 1);
end
