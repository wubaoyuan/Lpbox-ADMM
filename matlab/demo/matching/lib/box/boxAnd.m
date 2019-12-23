function box = boxAnd(boxs)
% Obtain the intersection of a set of boxes.
%
% Remark 1
%   The input "boxs" could be a cell array of boxes, e.g., boxs = {box1, box2, box3}, 
%   or a large matrix, e.g., boxs = [box1, box2, box3];
%
% Remark 2
%   Each box is a dim x 2 matrix, where box(i, 1) and box(i, 2) denote
%   the starting and ending position of the box in i-th dimension respectively.
%
% Remark 3
%   The size of the box could be calculated by 
%       siz = box(:, 2) - box(:, 1) + 1;
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

% cell -> double
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

% boundary of each dimension
box = zeros(dim, 2);
for di = 1 : dim
    box(di, 1) = max(boxs(di, 1, :));
    box(di, 2) = min(boxs(di, 2, :));
end
