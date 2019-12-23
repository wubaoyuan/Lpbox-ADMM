function dec = vis2dec(vis)
% Convert binary number to decimal.
%
% Example
%   in      -  vis = [0, 1, 0, 0]
%   call    -  dec = vis2dec(vis)
%   out     -  dec = 2 
%                  = 0 * 2 ^ 0 + 1 * 2 ^ 1 + 0 * 2 ^ 2 + 0 * 2 ^ 3
%
% Input
%   vis     -  binary vector, 1 x m
%
% Output
%   dec     -  decimal value
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-27-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

vis = vis(end : -1 : 1);
m = length(vis);
s = zeros(1, m);
for i = 1 : m
    s(i) = select(vis(i) == 1, '1', '0');
end
s = char(s);
dec = bin2dec(s);
