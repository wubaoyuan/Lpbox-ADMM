function Us = mkcca(Ks, par)
% Multi-set Kernel Canonical Correlation Analysis (mCCA).
%
% Remark
%   It can also be used to optimize KCCA for two sets of points.
%
% Input
%   Ks      -  original kernel matrix. 1 x m (cell), n x n
%   par     -  parameter
%     cen   -  flag of centrarization, {'y'} | 'n'
%     lams  -  manually specified regularization value, {0}
%     d     -  deciding the dimensions after projection, {.999}
%                d < 1:  percentage of energy to keep
%                1 <= d: #dimensions to keep
%
% Output
%   Us      -  cooeficient matrix, 1 x m (cell), n x d
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-20-2013

% function parameter
isCen = psY(par, 'cen', 'y');
lams = ps(par, 'lams', 0);
d = ps(par, 'd', .999);

% dimension
m = length(Ks);
ns = cellDim(Ks, 1);

% regularization
if isempty(lams)
    lams = zeros(1, m);
elseif length(lams) == 1
    lams = zeros(1, m) + lams;
end

% centralize
mes = cell(1, m);
for i = 1 : m
    if isCen
        [Xs{i}, mes{i}] = cenK(Ks{i});
    else
        mes{i} = zeros(ds(i), 1);
    end
end

% core implementation
Us = mkccaCore(Ks, ds, lams, d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Us = mkccaCore(Ks, ns, lams, d)
% Core implementation of MCCA.
%
% Input
%   Ks    -  original kernel matrix. 1 x m (cell), n x n
%   ns    -  dimension, 1 x m
%   lams  -  regularization weight, 1 x m
%   d     -  #dimension to keep
%
% Output
%   Us    -  cooeficient matrix, 1 x m (cell), n x d

% dimension
dSum = sum(ns);
m = length(ns);
s = n2s(ns);

% convariance
[CD, CO] = zeross(dSum, dSum);
for i = 1 : m
    % diagonal part
    idxi = s(i) : s(i + 1) - 1;
    CD(idxi, idxi) = Xs{i} * Xs{i}';

    % off-diagonal part
    for j = i + 1 : m
        idxj = s(j) : s(j + 1) - 1;
        CO(idxi, idxj) = Xs{i} * Xs{j}';
        CO(idxj, idxi) = CO(idxi, idxj)';
    end
end

% regularization
lam2s = lams;
for i = 1 : m
    idx = s(i) : s(i + 1) - 1;
    CD(idx, idx) = (1 - lam2s(i)) * CD(idx, idx) + lam2s(i) * eye(ds(i));
end

% generalized eigen-decomposition
[V, Lamb] = eig(CO, CD);
[lambs, idx] = sort(sum(Lamb, 2), 'descend');
V = V(:, idx);

% energy
lambs = lambs(1 : min(ds));
if d < 1
    d = thEgy(lambs, d);
end

% split
Vs = cell(1, m);
Us = cell(1, m);
for i = 1 : m
    Vs{i} = V(s(i) : s(i + 1) - 1, 1 : d);
    
    % constraint
    Us{i} = Vs{i}' * CD(s(i) : s(i + 1) - 1, s(i) : s(i + 1) - 1) * Vs{i};
end
U = sum(cat(3, Us{:}), 3);
