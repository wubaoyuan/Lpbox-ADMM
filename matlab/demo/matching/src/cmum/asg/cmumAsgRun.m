function wsRun = cmumAsgRun(tagSrc, tagAlg, varargin)
% Run graph matching algorithm on the CMU Motion data set.
%
% Input
%   tagSrc  -  source type, 1 | 2 | 3
%   tagAlg  -  algorithm type, 1 | 2 | ...
%   varargin
%     save option
%
% Output
%   wsRun
%     prex  -  name
%     Me    -  mean, nAlg x nBin
%     Dev   -  standard deviation, nAlg x nBin
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-04-2013

% save option
prex = cellStr('cmum', 'tagSrc', tagSrc, 'tagAlg', tagAlg);
[svL, path] = psSv(varargin, ...
                   'prex', prex, ...
                   'subx', 'run', ...
                   'fold', 'cmum/asg/run');

% load
if svL == 2 && exist(path, 'file')
    wsRun = matFld(path, 'wsRun');
    prInOut('cmumAsgRun', 'old, %s', prex);    
    return;
end
prIn('cmumAsgRun', 'new, %s', prex);

% parameters for generating src
[~, gaps, PFs] = cmumAsgPair(tagSrc);

% parameters for algorithms
[parAlgs, algs] = gmPar(tagAlg);

% dimension
nBin = length(gaps);
nReps = cellDim(PFs, 2);
nAlg = length(parAlgs);

% per gap
[Me, Dev, ObjMe, ObjDev] = zeross(nAlg, nBin);
prCIn('bin', nBin, 1);
for iBin = 1 : nBin
    prC(iBin);

    wsBin = cmumAsgRunBin(tagSrc, tagAlg, iBin, 'svL', 2);
    [Obj, Acc] = stFld(wsBin, 'Obj', 'Acc');
    Obj = Obj ./ repmat(Obj(end, :), nAlg, 1);
    
    % mean & deviation
    Me(:, iBin) = mean(Acc, 2);
    Dev(:, iBin) = std(Acc, 0, 2);
    ObjMe(:, iBin) = mean(Obj, 2);
    ObjDev(:, iBin) = std(Obj, 0, 2);
end
prCOut(nBin + 1);

% store
wsRun.prex = prex;
wsRun.Me = Me;
wsRun.Dev = Dev;
wsRun.ObjMe = ObjMe;
wsRun.ObjDev = ObjDev;

% save
if svL > 0
    save(path, 'wsRun');
end

prOut;
