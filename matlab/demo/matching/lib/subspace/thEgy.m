function [b, lamb2s] = thEgy(lambs, egy)
% Pick dimension by thresholding on the accumulated energy.
%   sum(lambs(1 : b)) / sum(lambs) >= egy
%
% Example
%   input   -  lambs = [1, .2, .1, .05]; 
%              egy = .95;
%   call    -  [b, lamb2s] = thEgy(lambs, egy)
%   output  -  b = 3;
%              lamb2s = [0.7407, 0.8889, 0.9630, 1];
%
% Input
%   lambs   -  energy vector (descendly sorted, positive), d x 1
%   egy     -  threshold
%
% Output
%   b       -  retained dimension
%   lamb2s  -  accumulated energy, d x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-12-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-16-2011

lambs = lambs / sum(lambs);
lamb2s = cumsum(lambs);
b = find(lamb2s >= egy, 1);

if isempty(b)
    b = length(lambs);
end
