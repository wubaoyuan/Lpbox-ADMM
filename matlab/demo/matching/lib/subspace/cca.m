function [v, Y1, Y2, V1, V2] = cca(X1, X2, par)
% Canonical Correlation Analysis (CCA).
%
% Input
%   X1      -  1st original sample matrix. (d1 + 1) x n
%   X2      -  2nd original sample matrix. (d2 + 1) x n
%   par     -  parameter
%     reg   -  method of regularization, 'fix' | {'cross'} | 'bach'
%              'fix': manually specified
%              'cross': cross-validation
%              'bach': lams is set to [n / 2; n / 2]
%                      See the paper in JMLR 2002, Kernel Independent Component Analysis, for details.
%     lams  -  manually specified regularization value (used if reg == 'fix'), {[0; 0]}
%     k     -  number of k in k-fold cross validation (used if reg == 'cross'), {10}
%     ms    -  number of grid search (used if reg == 'cross'), {[5; 5]}
%     mis   -  minimum value of lambda in search (used if reg == 'cross'), {[1e-4; 1e-4]}
%     mas   -  maximum value of lambda in search (used if reg == 'cross'), {[1e4; 1e4]}
%     sam   -  sampling method to create the grid (used if reg == 'cross'), {'log'} | 'lin'
%     egy   -  percentage of energy to keep for the choice of b, {[]}
%     b     -  #diemnsion kept after projection
%              if egy == [], then b has to be specified.
%     debg  -  debug flag, 'y' | {'n'}
%
% Output
%   v       -  objective value of CCA
%   Y1      -  1st transformed sample matrix, (b + 1) x n
%   Y2      -  2nd transformed sample matrix, (b + 1) x n
%   V1      -  1st transformation matrix, (d1 + 1) x (b + 1)
%   V2      -  2nd transformation matrix, (d2 + 1) x (b + 1)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function parameter
reg = ps(par, 'reg', 'bach');
lams = ps(par, 'lams', [0; 0]);
k = ps(par, 'k', 10);
ms = ps(par, 'ms', [10; 10]);
mis = ps(par, 'mis', [1; 1] * 1e-4);
mas = ps(par, 'mas', [1; 1] * 1e4);
sam = ps(par, 'sam', 'log');
egy = ps(par, 'egy', []);
b = ps(par, 'b', []);
isDebg = psY(par, 'debg', 'n');

% dimension
[d1, n] = size(X1); 
d2 = size(X2, 1); 

% using homographic coordinate in input
d1 = d1 - 1;
d2 = d2 - 1;
X10 = X1;
X20 = X2;
X1(end, :) = [];
X2(end, :) = [];

% centralize sample matrix
[X1, me1] = cenX(X1);
[X2, me2] = cenX(X2);

% convariance
C11 = X1 * X1';
C22 = X2 * X2';
C12 = X1 * X2';

% regularization
if rank(C11) < d1 || rank(C22) < d2
    if strcmp(reg, 'fix')
    elseif strcmp(reg, 'bach')
        lams = [1; 1] * n / 2;
    elseif strcmp(reg, 'cross')
        lams = ccaCross(X1, X2, k, ms, mis, mas, sam, egy, b, isDebg);
    else
        error('unknown regularization algorithm: %s', reg);
    end
    prom('t', 'reg %.5f %.5f\n', lams(1), lams(2));
end

% main algorithm of cca
[V1, V2] = ccaCore(C11, C22, C12, lams, egy, b);

% projection
Y1 = V1' * X1;
Y2 = V2' * X2;

% objective value
v = ccaObj(C11, C22, C12, V1, V2, lams);

% using homographic coordinate in output
Y1 = [Y1; ones(1, n)];
Y2 = [Y2; ones(1, n)];
V1 = [V1, zeros(d1, 1); - me1' * V1, 1];
V2 = [V2, zeros(d2, 1); - me2' * V2, 1];

% equal('Y1', Y1, V1' * X10);
% equal('Y2', Y2, V2' * X20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = ccaObj(C11, C22, C12, V1, V2, lams)
% Objective of CCA.
%
% Input
%   C11   -  covariance, d1 x d1
%   C22   -  covariance, d2 x d2
%   C12   -  covariance, d1 x d2
%   V1    -  1st transformation matrix, d1 x b
%   V2    -  2nd transformation matrix, d2 x b
%   lams  -  regularization weight, 2 x 1
%
% Output
%   Y1    -  1st transformed sample matrix, b x n
%   Y2    -  2nd transformed sample matrix, b x n
%   v     -  objective value of CCA

% global
global C11G C22G C12G;
isG = isempty(C11);

if isG
    [d1, d2] = size(C12G);
    A = V1' * C12G * V2;
    B1 = V1' * (C11G + lams(1) * eye(d1)) * V1;
    B2 = V2' * (C22G + lams(2) * eye(d2)) * V2;
    
    B1a = V1' * C11G * V1;
    B2a = V2' * C22G * V2;
else
    [d1, d2] = size(C12);
    A = V1' * C12 * V2;
    B1 = V1' * (C11 + lams(1) * eye(d1)) * V1;
    B2 = V2' * (C22 + lams(2) * eye(d2)) * V2;
    
    B1a = V1' * C11 * V1;
    B2a = V2' * C22 * V2;
end

% v = trace(A);

% v = trace((B1a + B2a) \ A);

a = trace(A);
b1 = trace(B1a);
b2 = trace(B2a);
if abs(b1) < 10e-10 || abs(b2) < 10e-10
    v = 0;
else
    v = a / (sqrt(b1) * sqrt(b2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V1, V2] = ccaCore(C11, C22, C12, lams, egy, b)
% Main implementation of CCA.
%
% Input
%   C11   -  covariance, d1 x d1
%   C22   -  covariance, d2 x d2
%   C12   -  covariance, d1 x d2
%   lams  -  regularization weight, 2 x 1
%   egy   -  energy threshold
%   b     -  #dimension to keep
%
% Output
%   V1    -  1st transformation matrix, d1 x b
%   V2    -  2nd transformation matrix, d2 x b

% global
global C11G C22G C12G;
isG = isempty(C11);

% covariance
if isG
    [d1, d2] = size(C12G);
    C1 = [zeros(d1, d1), C12G; ...
          C12G', zeros(d2, d2)];
    C2 = [C11G + lams(1) * eye(d1), zeros(d1, d2); ...
          zeros(d2, d1), C22G + lams(2) * eye(d2)];

else
    [d1, d2] = size(C12);
    C1 = [zeros(d1, d1), C12; ...
          C12', zeros(d2, d2)];
    C2 = [C11 + lams(1) * eye(d1), zeros(d1, d2); ...
          zeros(d2, d1), C22 + lams(2) * eye(d2)];
end

% generalized eigen-decomposition
% C1
% C2
[V, D] = eig(C1, C2);
[Lamb, index] = sort(sum(D, 2), 'descend');
V = V(:, index);

% [V2, D2] = eig(C1 + C2, C2);
% [Lamb2, index] = sort(sum(D2, 2), 'descend');
% V2 = V2(:, index);

% energy
d = min(d1, d2);
Lamb = Lamb(1 : d);
if ~isempty(egy)
    b = thEgy(Lamb, egy);
end

V1 = V(1 : d1, 1 : b);
V2 = V(d1 + 1 : end, 1 : b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lams = ccaCross(X1, X2, k, ms, mis, mas, sam, egy, b, isDebg)
% Cross validation for computing the regularization weights.
%
% Input
%   X1    -  1st original sample matrix. d1 x n
%   X2    -  2nd original sample matrix. d2 x n
%   k     -  the number of k in k-fold cross validation
%   ms    -  number of grid search if reg == 'cross', 2 x 1
%   mis   -  minimum value of lambda in search, 2 x 1
%   mas   -  maximum value of lambda in search, 2 x 1
%   sam   -  sampling method to create the grid, 'log' | 'lin'
%   egy   -  energy threshold
%   b     -  #dimension to keep
%
% Output
%   lams  -  regularization weights, 2 x 1

% global
global C11G C22G C12G;

% dimension
[~, n] = size(X1);

% randomly split data points into k fold
[~, idx] = divN(n, k);

% grid of regularization weight
grids = cell(1, 2);
for i = 1 : 2
    if strcmp(sam, 'log')
        mis(i) = log10(mis(i));
        mas(i) = log10(mas(i));
        grids{i} = logspace(mis(i), mas(i), ms(i));
    elseif strcmp(sam, 'lin')
        grids{i} = linspace(mis(i), mas(i), ms(i));
    else
        error(['unknown sampling method: ' sam]);
    end
end

% k-fold
[Objs, Bs] = zeross(ms(1), ms(2), k);
for c = 1 : k
    visTr = idx ~= c;
    visTe = idx == c;

    % sample for training
    X1Tr = X1(:, visTr);
    X2Tr = X2(:, visTr);
    X1Tr = cenX(X1Tr);
    X2Tr = cenX(X2Tr);

    % sample for testing
    X1Te = X1(:, visTe);
    X2Te = X2(:, visTe);
    X1Te = cenX(X1Te);
    X2Te = cenX(X2Te);
    
    % convariance
    C11G = X1Tr * X1Tr';
    C22G = X2Tr * X2Tr';
    C12G = X1Tr * X2Tr';
    C11GTe = X1Te * X1Te';
    C22GTe = X2Te * X2Te';
    C12GTe = X1Te * X2Te';

    % grid search for b
    for i1 = 1 : ms(1)
        for i2 = 1 : ms(2)
            
            % lambda
            lams = [grids{1}(i1); grids{2}(i2)];

            % cca
            V1 = ccaCore([], [], [], lams, egy, b);
            
            Bs(i1, i2, c) = size(V1, 2);
        end
    end
    
    % grid search for obj
    for i1 = 1 : ms(1)
        for i2 = 1 : ms(2)
            
            % lambda
            lams = [grids{1}(i1); grids{2}(i2)];

            % cca
            [V1, V2] = ccaCore([], [], [], lams, [], min(vec(Bs(:, :, c))));

            % objective value
            Objs(i1, i2, c) = ccaObj(C11GTe, C22GTe, C12GTe, V1, V2, lams);
        end
    end
end
Obj = sum(Objs, 3) / k;

% remove Nan
Obj(isnan(Obj)) = 0;

% debug
if isDebg
    rows = 1; cols = 2;
    axs = iniAx(10, rows, cols, [400 * cols, 400 * rows]);

    [GridX, GridY] = meshgrid(grids{1}, grids{2});
    
    set(gcf, 'CurrentAxes', axs{1});
    mesh(GridX, GridY, Obj);
    set(axs{1}, 'View', [135 45]);
    title('regularization');
    
    set(gcf, 'CurrentAxes', axs{2});
    contour(GridX, GridY, Obj);
    
    if strcmp(sam, 'log')
        set(axs{1}, 'XScale', 'log', 'YScale', 'log');
        set(axs{2}, 'XScale', 'log', 'YScale', 'log');
    else
        set(axs{1}, 'XScale', 'linear', 'YScale', 'linear');
        set(axs{2}, 'XScale', 'linear', 'YScale', 'linear');
    end
end

% optimum value
[tmp, p] = max(Obj(:));
[p1, p2] = ind2sub(ms, p);
lams = [grids{1}(p1); grids{2}(p2)];
