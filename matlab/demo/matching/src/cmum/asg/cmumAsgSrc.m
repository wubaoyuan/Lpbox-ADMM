function wsSrc = cmumAsgSrc(tag, pFs, nIns, varargin)
% Generate CMU Motion source for assignment problem.
%
% Input
%   tag      -  'house' | 'hotel'
%   pFs      -  frame index, 1 x 2
%   nIns     -  #inliers, 1 x 2, [1~30, 1~30]
%   varargin
%     save option
%
% Output
%   wsSrc
%     prex   -  prex
%     asgT   -  ground truth assignment
%     Pts    -  graph node set, 1 x mG (cell), 2 x ni
%     ords   -  order, 1 x mG (cell), 1 x ni
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 03-03-2013

% save option
prex = cellStr(tag, pFs, nIns);
[svL, path] = psSv(varargin, ...
                   'prex', prex, ...
                   'subx', 'src', ...
                   'fold', 'cmum/asg');

% load
if svL == 2 && exist(path, 'file')
    wsSrc = matFld(path, 'wsSrc');
    prInOut('cmumAsgSrc', 'old, %s', prex);
    return;
end
prIn('cmumAsgSrc', 'new, %s', prex);

% ground-truth label
CMUM = cmumHuman;

% marker position
Pts = CMUM.(tag).XTs(pFs);

% dimension
mG = 2;
n0 = size(Pts{1}, 2);
ns = nIns;

% per graph
ords = cellss(1, mG);
for iG = 1 : mG
    % randomly pick inlier nodes
    ord = randperm(n0);
    ord = ord(1 : nIns(iG));
    
    % store
    ords{iG} = ord;
    Pts{iG} = Pts{iG}(:, ord);
end

% ground-truth assignment
XT = eye(n0);
asgT.alg = 'truth';
asgT.X = XT(ords{1}, ords{2});

% store
wsSrc.prex = prex;
wsSrc.Pts = Pts;
wsSrc.asgT = asgT;
wsSrc.ords = ords;
wsSrc.tag = tag;
wsSrc.pFs = pFs;
wsSrc.nIns = nIns;

% save
if svL > 0
    save(path, 'wsSrc');
end

prOut;
