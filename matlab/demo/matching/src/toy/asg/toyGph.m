function Pt = toyGph(mes, Vars)
% Generate toy graph.
%
% Input
%   mes     -  mean of Gaussian, 1 x k (cell), 2 x 1
%   Vars    -  variance of Gaussian, 1 x k (cell), 2 x 2
%
% Output
%   Pt      -  graph node set, 2 x k
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-16-2012

if isempty(mes)
    Pt = [];
else
    Pt = smpGausses(mes, Vars, 1);
end
