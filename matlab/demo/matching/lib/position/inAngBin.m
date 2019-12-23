function iB = inAngBin(x, y, nB)
% Obtain the angle and put it into bins
%
% Input
%   x       -  x value, 1 x n
%   y       -  y value, 1 x n
%   nB      -  number of bins
%
% Output
%   iB      -  position of bins, 1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% angle
n = length(x);
P = atan2(y, x);

% starting angle for each bin
sp = linspace(-pi, pi, nB + 1);

% angle by angle
iB = zeros(1, n);
for i = 1 : n
    p = P(i);
    if p < sp(1)
        ib = 1;
    elseif p >= sp(nB + 1)
        ib = nB;
    else
        for j = 1 : nB
            if p >= sp(j) && p < sp(j + 1)
                ib = j;
                break;
            end
        end
    end
    
    iB(i) = ib;
end

