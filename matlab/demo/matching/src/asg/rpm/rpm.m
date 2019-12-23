function asg = rpm(gphs, asgT, par)
% Robust point matching (PRM).
%
% Reference
%   H. Chui and A. Rangarajan, "A new point matching algorithm for non-rigid registration",
%   Comput. Vis. Image Underst., 89(2-3):114-141, 2003.
%    
% Input
%   gphs     -  graphs, 1 x 2 (cell)
%     Pt     -  node position, 2 x ni
%   asgT     -  ground-truth assignment (can be [])
%   par      -  function parameter
%     algT   -  method for transformation, {'tps'} | 'rbf'
%     algX   -  method for computing correspondence, {'icp'} | 'mixture'
%     sigma  -  sigma, {1}
%     TFac   -  {500}
%     TRat   -  {.93}
%     deb    -  debug flag, 'y' | {'n'}
%
% Output
%   asg
%     alg    -  'rpm'
%     m      -  correspondence
%     vx     -  data points after warping, xmax x 2
%     acc    -  accuracy
%     obj    -  objective
%
% History
%   create   -  Anand Rangarajan (anand@noodle.med.yale.edu), 04-27-2000
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-03-2012

% function parameter
algT = ps(par, 'algT', 'tps');
algX = ps(par, 'algX', 'icp');
sigma = ps(par, 'sigma', 1);
TFac = ps(par, 'TFac', 500);
TRat = ps(par, 'TRat', 0.93);
isDeb = psY(par, 'deb', 'n');
prIn('rpm', 'algT %s, algX %s, sigma %.2f', algT, algX, sigma);

% init control parameters
lamda10 = 1;
lamda20 = 1;

% data points
x = gphs{2}.Pt';
y = gphs{1}.Pt';
vx = x;

% dimension
[xmax, dim] = size(x);
ymax = size(y, 1);

% outlier
T0 = max(x(:, 1)) ^ 2;
moutlier = 1 / sqrt(T0) * exp(-1);
m_outliers_row = ones(1, ymax) * moutlier;
m_outliers_col = ones(xmax, 1) * moutlier;

% init transformation parameters
theta = 0;
c_tps = zeros(xmax, dim + 1);
d_tps = eye(dim + 1, dim + 1);
w = zeros(xmax + dim + 1, dim);

% annealing procedure
T = T0;
TFin = T0 / TFac;
nItMa = ceil(log(TFin / T0) / log(TRat)) + 1;
nItMaT = 5;
if strcmp(algX, 'icp')
    nItMaT = 1;
end

pr('T0 %.4f, moutlier %.4f', T0, moutlier);

% debug
if isDeb
    rows = 2; cols = 3;
    Ax = iniAx(10, rows, cols, [300 * rows, 300 * cols]);
end

% optimization iteration
prCIn('anneal', nItMa, .1);
for nIt = 1 : nItMa
    prC(nIt);

    % repeat at each temperature
    for i = 1 : nItMaT
        % optimal correspondence
        X = rpmCorOpt(vx, y, T, algX, m_outliers_row, m_outliers_col, nIt, sigma);

        % optimal transformation
%        vy = X * y ./ ((sum(X'))' * ones(1, dim));
        vy = X * y ./ (X * ones(ymax, dim));
        %disp(min(X * ones(ymax, 1)));

        lamda1 = lamda10 * length(x) * T;
        lamda2 = lamda20 * length(x) * T;
%         lamda1 = lamda10;
%         lamda2 = lamda20;
        [c_tps, d_tps, w] = rpmTranOpt(algT, lamda1, lamda2, sigma, x, vy);

        % apply transformation
        vx = rpmWarp(algT, x, c_tps, d_tps, w, sigma);
    end

    % reduce temperature
    T = T * TRat;

    % stop condition
    if T < TFin
        break;
    end
    pr('T %.4f, lamda1 %.4f, lamda2 %.4f', T, lamda1, lamda2);

    % debug
    if isDeb
        rpmPlot2(Ax, x, y, vx, X, 1 / xmax, T, algT, c_tps, d_tps, w, sigma, algX);
        drawnow;
    end
end
prCOut(nIt + 1);

% continuous -> discrete
X = round(X);
X = X';

% matching with ground-truth assignment if possible
acc = matchAsg(X, asgT);

% transformation
tran = st('algT', algT, 'c_tps', c_tps, 'd_tps', d_tps, 'w', w);

% apply transformation
gphs{2}.Pt = tranRun(gphs{2}.Pt, tran);
gphs{1}.Eg = [];
gphs{2}.Eg = [];

% store
asg.alg = sprintf('%s_%s', algT, algX);
asg.X = X';
asg.acc = acc;
asg.tran = tran;
asg.gphs = gphs;
asg.obj = icpObj(gphs, asg.X, tran);

prOut;
