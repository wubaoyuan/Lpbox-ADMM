function wsFeat = cmumFeat(src, parF, varargin)
% Obtain the feature of CMU Motion sequence.
%
% Input
%   src     -  cmum src
%   parF    -  feature parameter
%   varargin 
%     save option
%
% Output
%   wsFeat
%     Pts   -  graph nodes, 1 x mG (cell), 2 x ni
%     Egs   -  graph edges, 1 x mG (cell), 2 x 2mi
%     Gs    -  node-edge adjacency, 1 x mG (cell), ni x mi
%     Hs    -  augumented node-edge adjacency, 1 x mG (cell), ni x (mi + ni)
%     XPs   -  node feature, 1 x mG (cell), dP x ni
%     XQs   -  edge feature, 1 x mG (cell), dQ x 2mi
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-31-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% save option
[svL, path] = psSv(varargin, 'fold', 'cmum/feat', ...
                             'prex', src.nm, ...
                             'subx', 'feat');

% load
if svL == 2 && exist(path, 'file')
    wsFeat = matFld(path, 'wsFeat');
    prInOut('cmumFeat', 'old, %s', src.nm);
    return;
end
prIn('cmumFeat', 'new, %s', src.nm);

% graph node
Pts = src.XTs;

% graph edge
[Egs, Gs, Hs] = gphConn(Pts, parF);

% graph feature
XPs = shpCtx(Pts, parF);
XQs = gphFeatQ(Pts, Egs);

% store
wsFeat.prex = src.nm;
wsFeat.Pts = Pts;
wsFeat.Egs = Egs;
wsFeat.Gs = Gs;
wsFeat.Hs = Hs;
wsFeat.XPs = XPs;
wsFeat.XQs = XQs;

% save
if svL > 0
    save(path, 'wsFeat');
end

prOut;
