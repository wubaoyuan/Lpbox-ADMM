function pathDStore(KP0, KQ0, Ct0, XQ10, XQ20, gphs, algX, XT0)
% Store the necessary variables for path-following algorithm.
%
% Input
%   KP0     -  node affinity matrix, n1 x n2
%   KQ0     -  edge affinity matrix, m1 x m2
%   Ct0     -  constraints, n1 x n2
%   XQ10    -  component 1, d x m1 | []
%   XQ20    -  component 2, d x m2 | []
%   gphs    -  graphs, 1 x 2 (cell)
%   algX    -  'ipfp' | 'path'
%   XT0     -  ground-truth correspondence, n1 x n2 | []
%
% History   
%   create  -  Feng Zhou (zhfe99@gmail.com), 09-01-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% global variable
global KP KQ Ct;
global G1 G2 H1 H2 GG1 GG2 HH1 HH2;
global G1s G2s H1s H2s HH1s HH2s;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T IndHH1 IndHH2;
global indG1 indG2 indH1 indH2;
global QQ1 QQ2 GHHQQG gamma;
global ns XT;

% store
KP = KP0;
KQ = KQ0;
Ct = Ct0;
XT = XT0;

% graph elements
[G1, H1] = stFld(gphs{1}, 'G', 'H');
[G2, H2] = stFld(gphs{2}, 'G', 'H');

% dimension
[n1, m1] = size(G1);
[n2, m2] = size(G2);
ns = [n1, n2];

% add additional nodes to make sure n1 == n2 (only for path-following)
if strcmp(algX, 'path')
    if n1 < n2
        mi = min(KP(:));
%        mi = 0;
        KP = [KP; zeros(n2 - n1, n2) + mi];
        G1 = [G1; zeros(n2 - n1, m1)];
        H1 = [H1; zeros(n2 - n1, m1)];
        Ct = [Ct; ones(n2 - n1, n2)];
        if ~isempty(XT)
            XT = [XT; zeros(n2 - n1, n2)];
        end
    elseif n1 > n2
        mi = min(KP(:));
%        mi = 0;
        KP = [KP, zeros(n1, n1 - n2) + mi];
        G2 = [G2; zeros(n1 - n2, m2)];
        H2 = [H2; zeros(n1 - n2, m2)];
        Ct = [Ct, ones(n1, n1 - n2)];
        if ~isempty(XT)
            XT = [XT, zeros(n1, n1 - n2)];
        end
    end
end

% % binary matrix -> index matrix (for saving memeory)
% [IndG1, IndG1T, IndH1, IndH1T] = mat2inds(G1, G1', H1, H1');
% [IndG2, IndG2T, IndH2, IndH2T] = mat2inds(G2, G2', H2, H2');

% auxiliary variables (for saving computational time)
GG1 = G1' * G1;
GG2 = G2' * G2;
HH1 = H1' * H1;
HH2 = H2' * H2;
[IndHH1, IndHH2] = mat2inds(HH1, HH2);

indG1 = mat2indC(G1);
indG2 = mat2indC(G2);
indH1 = mat2indC(H1);
indH2 = mat2indC(H2);

% sparse matrix
G1s = sparse(G1);
G2s = sparse(G2);
H1s = sparse(H1);
H2s = sparse(H2);
HH1s = H1s' * H1s;
HH2s = H2s' * H2s;

% factorize KQ using SVD
if isempty(XQ10)
    [U, S, V] = svd(KQ);
%     [U2, V2] = nmf_kl(KQ, 20);
%     max(max(abs(KQ - U2 * V2)))
    s = diag(S);
    idx = 1 : length(s);
    pr('svd: %d/%d', length(idx), length(s));
    k = length(idx);
    U = U(:, idx);
    V = V(:, idx);
    s = s(idx);

    XQ1 = multDiag('col', U, real(sqrt(s)));
    XQ2 = multDiag('col', V, real(sqrt(s)));
    XQ1 = XQ1';
    XQ2 = XQ2';
    
% already been factorized (eg, dgm)
else
    XQ1 = XQ10;
    XQ2 = XQ20;
end

% auxiliary variables for computing the derivative of the constant term
QQ1 = XQ1' * XQ1;
QQ2 = XQ2' * XQ2;
if strcmp(algX, 'path')
    GHHQQG = G1 * (HH1 .* QQ1) * G1' + G2 * (HH2 .* QQ2) * G2';
else
    GHHQQG = [];
end

gamma = multTr(QQ1 .* GG1 .* HH1) + multTr(QQ2 .* GG2 .* HH2);