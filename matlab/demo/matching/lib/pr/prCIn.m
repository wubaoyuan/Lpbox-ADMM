function prCIn(nm, nRep, sca)
% Start a propmter for counting.
%
% Input
%   nm      -  name
%   nRep    -  #steps
%   sca     -  scale of moving, (0, 1) | 1 | 2 | ...
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-29-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-03-2012

% variables set in "prSet.m"
global lPr;
global nmPrs ticPrs ticPr0s nRepPrs scaPrs;

% print
pr('%s: 0/%d', nm, nRep);

% insert
nmPrs{lPr} = nm;
ticPrs(lPr) = tic;
ticPr0s(lPr) = ticPrs(lPr);
nRepPrs(lPr) = nRep;

% scaling
if sca < 1
    sca = round(nRep * sca);
end
scaPrs(lPr) = sca;

lPr = lPr + 1;
