function wsFeat = cmumAsgFeat(wsSrc, parG, parF, varargin)
% Obtain the feature of CMU Motion for assignment problem.
%
% Input
%   wsSrc   -  cmum src
%   parG    -  graph parameter
%   parF    -  feature parameter
%   varargin
%     save option
%
% Output
%   wsFeat
%     gphs  -  graphs, 1 x mG (cell)
%     XPs   -  node feature, 1 x mG (cell), dP x ni
%     Fs    -  frames, 1 x mG (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-31-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

% save option
prex = wsSrc.prex;
[svL, path] = psSv(varargin, 'prex', prex, ...
                             'subx', 'feat', ...
                             'fold', 'cmum/asg');

% load
if svL == 2 && exist(path, 'file')
    wsFeat = matFld(path, 'wsFeat');
    prInOut('cmumAsgFeat', 'old, %s', prex);    
    return;
end
prIn('cmumAsgFeat', 'new, %s', prex);

% graph node
Pts = wsSrc.Pts;
mG = length(Pts);

% generate graphs
gphs = ps(wsSrc, 'gphs', []);
if isempty(gphs)
    gphs = newGphUs(Pts, parG);
end

% graph feature
XPs = shpCtxs(Pts, parF);

% frame
F0s = cmumVdo(wsSrc, 'svL', 2);
Fs = F0s(wsSrc.pFs);

% store
wsFeat.gphs = gphs;
wsFeat.XPs = XPs;
wsFeat.Fs = Fs;

% save
if svL > 0
    save(path, 'wsFeat');
end

prOut;
