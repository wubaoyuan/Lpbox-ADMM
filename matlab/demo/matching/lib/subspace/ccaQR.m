function [v, Y1, Y2, V1, V2, me1, me2] = ccaQR(X1, X2, C)
% Canonical Correlation Analysis (CCA). Implemented by QR + SVD.
%
% Remark:
%   if d > n, PCA is used to first reduce the dimension.
%
% Input
%   X1      -  orginal sample matrix, d1 x n1
%   X2      -  orginal sample matrix, d2 x n2
%   C       -  correspondance matrix (used iff n1 ~= n2), 2 x nC
%
% Output
%   v       -  correlation
%   Y1      -  transformed sample matrix, b x n1
%   Y2      -  transformed sample matrix, b x n2
%   V1      -  affine transformation matrix, d1 x b
%   V2      -  affine transformation matrix, d2 x b
%   me1     -  mean vectors, d1 x 1
%   me2     -  mean vectors, d2 x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% check number of samples
[d1, n1] = size(X1);
[d2, n2] = size(X2);
if ~exist('C', 'var')
    if n1 ~= n2
        error('A correspondance matrix is needed because of the unmatched number of samples');
    end
    C = [1 : n1; 1 : n2];
end
[W, W1, W2, D1, D2] = C2W(C);
dd1 = diag(D1);
dd2 = diag(D2);

% trival CCA
if n1 == 1 || n2 == 1
    v = 0;
    Y1 = 0; Y2 = 0;
    V1 = zeros(d1, 1); V2 = zeros(d2, 1);
    me1 = mean(X1, 2); me2 = mean(X2, 2);
    return;
end

ns = [n1, n2];
ds = [d1, d2];
X0s = {X1, X2};
dds = {dd1, dd2};
[Xs, Qs, Rs, perms, mes, V2s] = cellss(1, 2);
rs = zeros(1, 2);
for i = 1 : 2
    n = ns(i);
    d = ds(i);
    X0 = X0s{i};
    dd = dds{i};

    % remove unused samples
    idx = find(dd == 0);
    if ~isempty(idx)
        X0(:, idx) = [];
        C(i, :) = C(i, :) - C(i, 1) + 1;
        n = n - length(idx);
        dd(idx) = [];
        dds{i} = dd;
    end

    % zero mean
    me = X0 * dd / sum(dd);
    X = X0 - repmat(me, 1, n);
    mes{i} = me;
    
    % reduce dimension by PCA
    if d > n
        Cov = (X .* repmat(dd', d, 1)) * X';
        [V, sig] = eigk(Cov, n - 1);
        r = sum(abs(sig) > eps);

        if r == 0
            v = 0;
            Y1 = 0; Y2 = 0;
            V1 = zeros(d1, 1); V2 = zeros(d2, 1);
            me1 = mean(X1, 2); me2 = mean(X2, 2);
            return;
        end

        V = V(:, 1 : r);
        X = (V' * X);
        d = r;
        V2s{i} = V;
    end

    % QR
    dd2 = dd .^ (.5);
    [Q0, R0, perms{i}] = qr(repmat(dd2, 1, d) .* X', 0);
    r = sum(abs(diag(R0)) > 10e-10); rs(i) = r;
    Qs{i} = Q0(:, 1 : r); Rs{i} = R0(1 : r, 1 : r);
    Xs{i} = X;
end
ra = min(rs);

% svd A = Q1' * (DH1 \ W / DH2) * Q2;
Q1 = repmat(dds{1} .^ (-.5), 1, rs(1)) .* Qs{1};
Q2 = repmat(dds{2} .^ (-.5), 1, rs(2)) .* Qs{2};
A = Q1(C(1, :), :)' * Q2(C(2, :), :);
[U1, Sig, U2] = svd(A, 0);
Us = {U1(:, 1 : ra), U2(:, 1 : ra)};

% correlation
v = diag(Sig(1 : ra, 1 : ra));
if nargout == 1
    return;
end

% restore
[Vs, Ys] = cellss(1, 2);
for i = 1 : 2
    perm = perms{i};
    r = rs(i);
    d = ds(i);
    X = Xs{i};
    
    V = Rs{i} \ Us{i};

    V(perm, :) = [V; zeros(d - r, ra)];
    
    Y = V' * X;

    Vs{i} = V;
    Ys{i} = Y;
    
    % PCA
    if ~isempty(V2s{i})
        Vs{i} = V2s{i} * V;
    end
end

Y1 = Ys{1}; Y2 = Ys{2};
V1 = Vs{1}; V2 = Vs{2};
me1 = mes{1}; me2 = mes{2};
