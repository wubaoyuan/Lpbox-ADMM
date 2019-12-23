function [mks, cls] = genMkCl(c)
% Return with markers and colors.
%
% Input
%   c       -  index (optional)
%
% Output
%   mks     -  string | 1 x 12 (cell) (if c has not been specified)
%   cls     -  string | 1 x 12 (cell) (if c has not been specified)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-03-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

mks = {'o', 's', '^', 'd', '+', '*', 'x', 'p', 'v', 'o', 's', '^', 'd'};
cls = {[1 0 0], [0 0 1], [0 1 0], [1 0 1], [0 0 0], [0 1 1], [.3 .3 .3], [.5 .5 .5], [.7 .7 .7], [.1 .1 .1], [1 .8 0], [1, .4, .6]};
%          'r',     'b',     'g',     'm',     'k',     'c',  

if exist('c', 'var')
    cc = mod(c - 1, length(cls)) + 1;
    mks = mks{cc};
    cls = cls{cc};
end
