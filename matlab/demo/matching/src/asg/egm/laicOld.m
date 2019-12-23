function asg = laic(Rs, gphTr, PtTe, Ct, XT, par)
% LAIC algorithm.
%
% Input
%   Rs        -  response, h x w x nTr
%   gphTr     -  training graph
%   PtTe      -  testing point set, 2 x nTe
%   EgTr      -  edge matrix for trainint point set, 2 x mTr
%   Ct        -  constraint, nTr x nTe
%   XT        -  ground-truth correspondence, nTr x nTe | []
%   par       -  parameter
%
% Output
%   asg
%     alg     -  'laic'
%     X       -  correspondence Hmatrix, Nm x Nt
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 12-13-2010
%   modify    -  Feng Zhou (zhfe99@gmail.com), 07-26-2012

% dimension
[nTr, nTe] = size(Ct);

% training graph
[PtTr, EgTr] = stFld(gphTr, 'Pt', 'Eg');

% edge -> adjacency matrix
A0 = gphEg2Adj(EgTr, size(PtTr, 2));

% make sure every node has three neighbours
valTr = laicValE(A0);
valTe = sum(Ct(valTr, :)) > 0;
if any(~valTr)
    pr('some point missed in graph');
end

% node similarity
K = zeros(nTr, nTe);
for iTr = 1 : nTr
    K(iTr, :) = log(imgPtVal(Rs(:, :, iTr), PtTe)');
end

% with constraints
K = K .* Ct;

% normalize
mi = min(K(Ct == 1));
ma = max(K(Ct == 1));
%K = (K - mi) / (ma - mi);

% simiarlity -> distance
D0 = -K;
D0(Ct == 0) = 1000;

% adjust the matrix
A = A0(valTr, valTr);
D = D0(valTr, valTe);
Pts = {PtTr(:, valTr), PtTe(:, valTe)};

% run
[X, Xs, C, H, idxTrs, idxTes] = laicCore(Pts, A, D, par);

% objective
obj = laicObj(X, Pts, C, H);
if ~isempty(XT)
    [objT, objT1, objT2] = laicObj(XT, Pts, C, H);
else
    objT = 0;
    objT1 = 0;
    objT2 = 0;
end
mX = length(Xs);
[objs, obj1s, obj2s] = zeross(1, mX);
for iX = 1 : mX
    [objs(iX), obj1s(iX), obj2s(iX)] = laicObj(Xs{iX}, Pts, C, H);
end

% restore to the original result
X0 = zeros(nTr, nTe);
X0(valTr, valTe) = X;
idxTrInv = find(valTr == 0);
for i = 1 : length(idxTrInv)
    iTr = idxTrInv(i);
    
    idx = find(Ct(iTr, :) == 1);
    [~, p] = min(D0(iTr, idx));
    iTe = idx(p);
    
    X0(iTr, iTe) = 1;
    idxTrs = [idxTrs; iTr];
    idxTes = [idxTes, iTe];
end

% store
asg.alg = 'laic';
asg.X = X0;
asg.Xs = Xs;
asg.idxTrs = idxTrs;
asg.idxTes = idxTes;
asg.obj = obj;
asg.objT = objT;
asg.objT1 = objT1;
asg.objT2 = objT2;
asg.objs = objs;
asg.obj1s = obj1s;
asg.obj2s = obj2s;