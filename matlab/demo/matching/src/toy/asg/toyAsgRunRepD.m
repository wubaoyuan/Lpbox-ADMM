function wsRep = toyAsgRunRepD(tagSrc, tagAlg, iRep, varargin)
% Run graph matching algorithm on the toy data set.
%
% Remark
%   Edge feature is asymmetric.
%
% Input
%   tagSrc   -  source type, 1 | 2 | 3
%   tagAlg   -  algorithm type, 1 | 2 | ...
%   iRep     -  #repetition index
%   varargin
%     save option
%
% Output
%   wsRep
%     prex   -  name
%     Acc    -  accuracy, nAlg x nBin
%     Obj    -  objective, nAlg x nBin
%            
% History    
%   create   -  Feng Zhou (zhfe99@gmail.com), 01-25-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% save option
prex = cellStr('toy', 'tagSrc', tagSrc, 'tagAlg', tagAlg, 'iRep', iRep);
[svL, path] = psSv(varargin, ...
                   'prex', prex, ...
                   'subx', 'rep', ...
                   'fold', 'toy/asg/repD');

% load
if svL == 2 && exist(path, 'file')
    wsRep = matFld(path, 'wsRep');
    prInOut('toyAsgRunRepD', 'old, %s', prex);
    return;
end
prIn('toyAsgRunRepD', 'new, %s', prex);

% parameters for generating src
ParSrc = toyAsgPar(tagSrc);
parSrcs = ParSrc(iRep, :);

% parameters for algorithms
[parAlgs, algs] = gmPar(tagAlg);

% dimension
nBin = length(parSrcs);
nAlg = length(parAlgs);

% per bin (pair)
[Acc, Obj] = zeross(nAlg, nBin);
prCIn('bin', nBin, 1);
for iBin = 1 : nBin
    prC(iBin);
    
    % src
    wsSrc = toyAsgSrcD(parSrcs{iBin}{:});
    [gphs, asgT] = stFld(wsSrc, 'gphs', 'asgT');

    % affinity
    parKnl = st('nor', 'n', 'alg', 'toy');
    [KP, KQ] = conKnlGphPQD(gphs, parKnl);
    K = conKnlGphKD(KP, KQ, gphs);
    Ct = ones(size(KP));

    % directed -> undirected
    gphUs = gphD2Us(gphs);
    [KU, KQU] = knlGphKD2U(KP, KQ, gphUs);

    % per algorithm
    for iAlg = 1 : nAlg
        % parameter
        pars = parAlgs{iAlg};

        if strcmpi(algs{iAlg}, 'PM')
            asg = pm(K, KQU, gphUs, asgT);
        elseif strcmpi(algs{iAlg}, 'fgm') || strcmpi(algs{iAlg}, 'fgm-u')
            asg = fgmU(KP, KQU, Ct, gphUs, asgT, pars{:});
            X = asg.X;
            asg.obj = X(:)' * K * X(:);
        elseif strcmpi(algs{iAlg}, 'fgm-d')
            asg = fgmD(KP, KQ, Ct, gphs, asgT, pars{:});
        else
            asg = gm(K, Ct, asgT, pars{:});
        end

        % objective
        Acc(iAlg, iBin) = asg.acc;
        Obj(iAlg, iBin) = asg.obj;
    end
end
prCOut(nBin + 1);

% store
wsRep.prex = prex;
wsRep.Acc = Acc;
wsRep.Obj = Obj;

% save
if svL > 0
    save(path, 'wsRep');
end

prOut;
