function wsSrc = toyAsgTraSrc(tag, k0, kNs, nF, headMa, varargin)
% Generate source of toy trajectory for assignment problem.
%
% Input
%   tag       -  shape of latent graph, 1 | ...
%                1: random graph
%   k0        -  #nodes
%   kNs       -  #noisy node, 1 x m
%   nF        -  #frames
%   headMa    -  maximum head position, 1 | 2 | ... | nF
%   varargin
%     save option
%
% Output
%   wsSrc
%     prex    -  prex
%     asgT    -  ground truth assignment
%     messs   -  mean of Gaussian, 1 x m (cell), 1 x nF (cell), 1 x ki (cell), d x 1
%     Varsss  -  variance of Gaussian, 1 x m (cell), 1 x nF (cell), 1 x ki (cell), d x 2
%     As      -  trajectory status, 1 x m (cell), ki x nF
%     Trass   -  trajectory position, 1 x m (cell), 1 x ki (cell), d x nF
%     Ptss    -  graph node set, 1 x m (cell), 1 x nF (cell), d x k_{t_i}^i
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% save option
prex = cellStr(tag, k0, kNs, nF, headMa);
[svL, path] = psSv(varargin, 'prex', prex, ...
                             'subx', 'src', ...
                             'fold', 'toy/asgTra');

% load
if svL == 2 && exist(path, 'file')
    prIn('toyAsgTraSrc', 'old, %s', prex);
    wsSrc = matFld(path, 'wsSrc');
    prOut;
    return;
end
prIn('toyAsgTraSrc', 'new, %s', prex);

% dimension
d = 2;
m = 2;
ks = k0 + kNs;

% model for node
[me0ss, Var0ss] = toyGphTraMod(tag, k0, nF);

% per graph
[messs, Varsss, As, Trass, Ptss, ords] = cellss(1, m);
for i = 1 : m
    % model for noisy node
    [meNss, VarNss] = toyGphTraMod(1, kNs(i), nF);
    
    % re-order
    ords{i} = randperm(ks(i));
    
    % per frame
    ki = ks(i);
    [messs{i}, Varsss{i}] = cellss(1, nF);
    X = zeros(d, ki, nF);
    for iF = 1 : nF
        % combine
        mes = cellCat(me0ss{iF}, meNss{iF});
        Vars = cellCat(Var0ss{iF}, VarNss{iF});
        
        % re-order
        messs{i}{iF} = mes(ords{i});
        Varsss{i}{iF} = Vars(ords{i});
        
        % generate node
        X(:, :, iF) = toyGph(messs{i}{iF}, Varsss{i}{iF});
    end

    % missing node in time
    As{i} = zeros(ki, nF);
    for c = 1 : ki
        % random position for head
        head = randi([1 headMa]);

        % random position for tail
        tail = nF;

        As{i}(c, head : tail) = 1;
    end

    % re-convert
    Trass{i} = mcvXA2Tra(X, As{i});
    Ptss{i} = mcvXA2Pt(X, As{i});
end

% ground-truth assignment
asgT.alg = 'truth';
PT = zeros(ks);
idx = sub2ind(ks, 1 : k0, 1 : k0);
PT(idx) = 1;
asgT.P = PT(ords{1}, ords{2});

% store
wsSrc.prex = prex;
wsSrc.asgT = asgT;
wsSrc.messs = messs;
wsSrc.Varsss = Varsss;
wsSrc.As = As;
wsSrc.Trass = Trass;
wsSrc.Ptss = Ptss;
wsSrc.ords = ords;

% save
if svL > 0
    save(path, 'wsSrc');
end

prOut;
