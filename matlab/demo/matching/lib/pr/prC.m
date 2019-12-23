function prC(iRep)
% Prompt information of a counter.
%
% Input
%   iRep    -  current step
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-29-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-02-2012

% variables set in "prSet.m"
global lPr;
global nmPrs ticPrs nRepPrs scaPrs;

lPr = lPr - 1;
if (iRep - 1 ~= 0 && mod(iRep - 1, scaPrs(lPr)) == 0) || (iRep - 1 == nRepPrs(lPr))
    % time
    t = toc(ticPrs(lPr));

    % print
    pr('%s: %d/%d, %.2f secs', nmPrs{lPr}, iRep - 1, nRepPrs(lPr), t);

    % re-start a timer
    ticPrs(lPr) = tic;
end    
lPr = lPr + 1;