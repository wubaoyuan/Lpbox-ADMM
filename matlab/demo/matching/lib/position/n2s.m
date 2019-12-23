function [s, H] = n2s(ns)
% Obtain segment information from the segment length.
%
% Example
%   input   -  ns = [4 2 3];
%   call    -  [s, H] = n2s(ns);
%   output  -  s = [1 5 7 10];
%              H = [1 1 1 1 0 0 0 0 0; ...
%                   0 0 0 0 1 1 0 0 0; ...
%                   0 0 0 0 0 0 1 1 1];
%
% Input
%   ns      -  segment length, 1 x m
%
% Output
%   s       -  starting position of each segment, 1 x (m + 1)
%   H       -  frame-segment indicator matrix, m x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-30-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

s = cumsum([1 ns(:)']);

if nargout == 2
    m = length(ns);
    n = sum(ns);
    H = zeros(m, n);
    for i = 1 : m
        H(i, s(i) : s(i + 1) - 1) = 1;
    end
end
