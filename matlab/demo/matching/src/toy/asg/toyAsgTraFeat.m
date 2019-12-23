function wsFeat = toyAsgTraFeat(wsSrc, parF, varargin)
% Obtain feature of toy trajectory for assignment problem.
%
% Input
%   wsSrc    -  toy src
%   parF     -  feature parameter
%   varargin
%     save option
%
% Output
%   wsFeat
%     Ptss   -  graph nodes, 1 x m (cell), 1 x nF (cell), 2 x ki
%     Egss   -  graph edges, 1 x m (cell), 1 x nF (cell), 2 x li
%     XPtss  -  node feature, 1 x m (cell), 1 x nF (cell), dPt x ki
%     XEgss  -  edge feature, 1 x m (cell), 1 x nF (cell), dEg x li 
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 07-31-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% save option
prex = wsSrc.prex;
[svL, path] = psSv(varargin, 'prex', prex, ...
                             'subx', 'feat', ...
                             'fold', 'toy/asg');

% load
if svL == 2 && exist(path, 'file')
    prIn('toyAsgFeat', 'old, %s', prex);
    wsFeat = matFld(path, 'wsFeat');
    prOut;
    return;
end
prIn('toyAsgFeat', 'new, %s', prex);

% graph node
Ptss = wsSrc.Ptss;

% dimension
m = length(Ptss);

% per graph
[Egss, XPtss, XEgss] = cellss(1, m);
for i = 1 : m
    % graph edge
    Egss{i} = gphPt2Eg(Ptss{i}, parF);
    
    % graph feature
    [XPtss{i}, XEgss{i}] = gphFeat(Ptss{i}, Egss{i}, parF);
end

% store
wsFeat.prex = prex;
wsFeat.Ptss = Ptss;
wsFeat.Egss = Egss;
wsFeat.XPtss = XPtss;
wsFeat.XEgss = XEgss;

% save
if svL > 0
    save(path, 'wsFeat');
end

prOut;
