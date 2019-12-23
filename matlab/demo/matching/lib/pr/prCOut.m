function t = prCOut(nRep)
% Stop a propmter for counting.
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-29-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-02-2012

% variables set in "prSet.m"
global lPr;
global nmPrs ticPrs ticPr0s nRepPrs scaPrs;

lPr = lPr - 1;

% time
t = toc(ticPr0s(lPr));

% print
pr('%s: %d iters, %.2f secs', nmPrs{lPr}, nRep, t);
