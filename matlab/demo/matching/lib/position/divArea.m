function Areas = divArea(box, rows, cols)
% Divide an box area into blocks.
%
% Input
%   box     -  box position, 1 x 4
%   rows    -  number of rows
%   cols    -  number of rows
%
% Output
%   Areas   -  the position of areas, 4 x rows x cols
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-14-2010

ave = floor(siz ./ [rows, cols]);
Areas = zeros(4, rows, cols);

for r = 1 : rows
    r1 = ave(1) * (r - 1) + 1;
    r2 = ave(1) * r;
    if r == rows
        r2 = siz(1); 
    end
    
    for c = 1 : cols
        c1 = ave(2) * (c - 1) + 1;
        c2 = ave(2) * c;
        if c == cols
            c2 = siz(2);
        end

        Areas(:, r, c) = [r1, c1, r2 - r1 + 1, c2 - c1 + 1]';
    end
end

