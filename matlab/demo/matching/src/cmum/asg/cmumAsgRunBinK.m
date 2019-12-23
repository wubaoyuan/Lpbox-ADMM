function wsK = cmumAsgRunBinK(tagSrc, iBin, varargin)
% Run graph matching algorithm on the CMUM Motion data set.
%
% Input
%   tagSrc  -  source type, 1 | 2 | 3
%   iBin    -  bin index
%   varargin
%     save option
%
% Output
%   wsK
%     Ks    -  affinity K, 1 x nRep (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-05-2013

% save option
prex = cellStr('cmum', 'tagSrc', tagSrc, 'iBin', iBin);
[svL, path] = psSv(varargin, ...
                   'prex', prex, ...
                   'subx', 'binK', ...
                   'fold', 'cmum/asg/bin');

% load
if svL == 2 && exist(path, 'file')
    wsK = matFld(path, 'wsK');
    prInOut('cmumAsgRunBinK', 'old, %s', prex);    
    return;
end
prIn('cmumAsgRunBinK', 'new, %s', prex);

% parameters for generating src
[tag, gaps, PFs, nIns] = cmumAsgPair(tagSrc);
gap = gaps(iBin);
PF = PFs{iBin};

% dimension
nRep = size(PF, 2);

% per repetition (pair)
[Ks, KPs, KQs, eis] = cellss(1, nRep);
[n1s, n2s, m1s, m2s, ranKs, ranQs] = zeross(1, nRep);
prCIn('nRep', nRep, 1);
for iRep = 1 : nRep
    prC(iRep);

    % src
    pFs = PF(:, iRep);
    wsSrc = cmumAsgSrc(tag, pFs, nIns);
    asgT = wsSrc.asgT;

    % feature
    parG = st('link', 'del');
    parF = st('smp', 'n', 'nBinT', 4, 'nBinR', 3);
    wsFeat = cmumAsgFeat(wsSrc, parG, parF);
    [gphs, XPs, Fs] = stFld(wsFeat, 'gphs', 'XPs', 'Fs');

    % affinity
    parKnl = st('alg', 'cmum');
    [KP, KQ] = conKnlGphPQS(gphs, parKnl);
    K = conKnlGphKS(KP, KQ, gphs);
    Ct = ones(size(KP));

    % symmetric graph -> asymmetric graph
    gphAs = gphS2As(gphs);
    KQA = [KQ, KQ; KQ, KQ];
    
    % rank & eigen values
    K = full(K);
    K = K - diag(diag(K));

    % store
    [n1s(iRep), n2s(iRep)] = size(KP);
    [m1s(iRep), m2s(iRep)] = size(KQA);
    ranKs(iRep) = rank(K);
    ranQs(iRep) = rank(KQA);
    eis{iRep} = eig(K);
end
prCOut(nRep + 1);

% store
wsK.prex = prex;

wsK.n1s = n1s;
wsK.n2s = n2s;
wsK.m1s = m1s;
wsK.m2s = m2s;
wsK.ranKs = ranKs;
wsK.ranQs = ranQs;
wsK.eis = eis;

% save
if svL > 0
    save(path, 'wsK');
end

prOut;
