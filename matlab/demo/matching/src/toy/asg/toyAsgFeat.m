function wsFeat = toyAsgFeat(wsSrc, parF, varargin)
% Obtain the feature of toy graph for assignment.
%
% Input
%   wsSrc   -  toy src
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
prex = wsSrc.prex;
[svL, path] = psSv(varargin, 'prex', prex, ...
                             'subx', 'feat', ...
                             'fold', 'toy/asg');

% load
if svL == 2 && exist(path, 'file')
    wsFeat = matFld(path, 'wsFeat');
    prInOut('toyAsgFeat', 'old, %s', prex);    
    return;
end
prIn('toyAsgFeat', 'new, %s', prex);

% graph node
Pts = wsSrc.Pts;

% graph edge
Egs = ps(wsSrc, 'Egs', []);
if isempty(Egs)
    [Egs, Gs, Hs] = gphConn(Pts, parF);
end

% graph feature
XPs = shpCtx(Pts, parF);
XQs = gphFeatQ(Pts, Egs);

% store
wsFeat.prex = prex;
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
