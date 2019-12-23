function [Ctxs, Pts] = shpCtxs(Pt0s, par)
% Compute shape context for point cloud.
%
% Input
%   Pt0s     -  original sequence of point set, 1 x nF (cell), 2 x nPt0i
%   par      -  function parameter, see shpCtx.m for more details
%
% Output
%   Ctxs     -  shape context, 1 x nF (cell), (nBinT x nBinR) x nPti
%   Pts      -  new sequence of point sets, 1 x nF (cell), 2 x nPti
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 09-26-2010
%   modify   -  Feng Zhou (zhfe99@gmail.com), 01-29-2012

% dimension
nF = length(Pt0s);
prIn('shpCtxs', 'nF %d', nF);

% per graph
prCIn('gph', nF, .1);
[Pts, Ctxs] = cellss(1, nF);
for iF = 1 : nF
    prC(iF);
    [Ctxs{iF}, Pts{iF}] = shpCtx(Pt0s{iF}, par);
end
prCOut(nF);

prOut;
