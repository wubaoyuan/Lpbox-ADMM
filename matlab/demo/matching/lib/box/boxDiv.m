function boxs = boxDiv(box0, ms)
% Divide the box into blocks.
%   Each box is a dim x 2 matrix, where 1st and 2nd column denote the starting and ending position respectively.
%
% Input
%   box0    -  box, dim x 2
%   ms      -  number of bins, dim x 1
%
% Output
%   boxs    -  box set, dim x 2 x ms(1) x ms(2) x ...
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-22-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

dim = size(box0, 1);
lens = floor((box0(:, 2) - box0(:, 1)) ./ ms);

boxs = zeros(dim, 2, ms(1), ms(2));

[head, tail] = zeross(2, 1);
for r = 1 : ms(1)
    head(1) = lens(1) * (r - 1) + box0(1, 1);
    tail(1) = select(r < ms(1), head(1) + lens(1), box0(1, 2));

    for c = 1 : ms(2)
        head(2) = lens(2) * (c - 1) + box0(2, 1);
        tail(2) = select(c < ms(2), head(2) + lens(2), box0(2, 2));

        boxs(:, :, r, c) = [head, tail];
    end
end
