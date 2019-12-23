function prSet(lMa)
% Set the promption level.
%
% Input
%   lMa     -  maximum level, 0 | 1 | 2 | ...
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-29-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-02-2012

% variables
global lPr lMaPr;
global nmPrs ticPrs ticPr0s nRepPrs scaPrs;

% level
lPr = 1;
lMaPr = lMa;

% list
nMa = 10;
nmPrs = cell(1, nMa);
ticPrs = zeros(1, nMa, 'uint64');
ticPr0s = zeros(1, nMa, 'uint64');
[nRepPrs, scaPrs] = zeross(1, nMa);
