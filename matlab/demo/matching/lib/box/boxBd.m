function box = boxBd(boxs)
% Obtain the bounding box of a set of boxes.
% 
% Remark
%   Each box is a dim x 2 matrix, where 1st and 2nd column denote
%   the minimum and maximum value of each dimension respectively.
%
% Input
%   boxs    -  set of bounding boxs, 1 x m (cell) or dim x 2 x m (double)
%
% Output
%   box     -  new bounding box, dim x 2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if iscell(boxs)
    m = length(boxs);
    dim = size(boxs{1}, 1);
    box0s = boxs;
    
    boxs = zeros(dim, 2, m);
    for i = 1 : m
        boxs(:, :, i) = box0s{i};
    end
else
    dim = size(boxs, 1);
end

box = zeros(dim, 2);
for di = 1 : dim
    box(di, 1) = min(boxs(di, 1, :));
    box(di, 2) = max(boxs(di, 2, :));
end
