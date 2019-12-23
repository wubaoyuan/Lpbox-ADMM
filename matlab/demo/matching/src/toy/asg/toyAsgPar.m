function [Par, vals] = toyAsgPar(tag)
% Obtain parameters for generating toy sources used for testing assignment algorithm.
%
% Input
%   tag     -  type of pair, 1 | 2 | 3
%                1 : different #outlier
%                2 : different deformation
%                3 : different density
%
% Output
%   Par     -  parameters, nRep x nBin (cell)
%   vals    -  value that are changing, 1 x nBin
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-07-2013

% different #outliers
if tag == 1

    nIn = 20;
    nOuts = 0 : 2 : 20;
    egDen = 1;
    egDef = 0;
    nBin = length(nOuts);
    nRep = 100;
    Par = cell(nRep, nBin);

    for iBin = 1 : nBin
        nOut = nOuts(iBin);
        for iRep = 1 : nRep
            Par{iRep, iBin} = {1, nIn, [0 0] + nOut, egDen, egDef};
        end
    end
    vals = nOuts;
    
% different #outliers
elseif tag == 4

    nIn = 20;
    nOuts = 0 : 2 : 20;
    egDen = 1;
    egDef = 0;
    nBin = length(nOuts);
    nRep = 20;
    Par = cell(nRep, nBin);

    for iBin = 1 : nBin
        nOut = nOuts(iBin);
        for iRep = 1 : nRep
            Par{iRep, iBin} = {1, nIn, [0 0] + nOut, egDen, egDef};
        end
    end
    vals = nOuts;

% different deformation
elseif tag == 2

    nIn = 20;
    nOut = 0;
    egDen = 1;
    egDefs = 0 : .02 : .2;
    nBin = length(egDefs);
    nRep = 100;
    Par = cell(nRep, nBin);

    for iBin = 1 : nBin
        egDef = egDefs(iBin);
        for iRep = 1 : nRep
            Par{iRep, iBin} = {1, nIn, [0 0] + nOut, egDen, egDef};
        end
    end
    vals = egDefs;

% different density
elseif tag == 3

    nIn = 20;
    nOut = 0;
    egDens = 1 : -.1 : .3;
    egDef = 0;
    nBin = length(egDens);
    nRep = 100;
    Par = cell(nRep, nBin);

    for iBin = 1 : nBin
        egDen = egDens(iBin);
        for iRep = 1 : nRep
            Par{iRep, iBin} = {1, nIn, [0 0] + nOut, egDen, egDef};
        end
    end
    vals = egDens;

else
    error('unknown tag: %d', tag);
end
