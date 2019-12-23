function ipfpUStore(KP0, KQ0, Ct0, gphs)
% Store the necessary variables for IPFP algorithm.
%
% Input
%   KP0     -  node affinity matrix, n1 x n2
%   KQ0     -  edge affinity matrix, m1 x m2
%   Ct0     -  constraints, n1 x n2
%   gphs    -  graphs, 1 x 2 (cell)
%     GA    -  node-edge adjacency, n x mi
%     HA    -  node-edge adjacency, n x (mi + ni)
%
% History   
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-14-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% global variable
global L KP KQ Ct;
global G1 G2 H1 H2 GG1 GG2 HH1 HH2;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T;

KP = KP0;
KQ = KQ0;
Ct = Ct0;

% graph elements
[G1, H1] = stFld(gphs{1}, 'G', 'H');
[G2, H2] = stFld(gphs{2}, 'G', 'H');

% dimension
[n1, m1] = size(G1);
[n2, m2] = size(G2);

% make sure n1 == n2
% if n1 < n2
%     KP = [KP; zeros(n2 - n1, n2)];
%     G1 = [G1; zeros(n2 - n1, m1)];
%     H1 = [G1, eye(n2)];
%     Ct = [Ct; ones(n2 - n1, n2)];
% elseif n1 > n2
%     KP = [KP, zeros(n1, n1 - n2)];
%     G2 = [G2; zeros(n1 - n2, m2)];
%     H2 = [G2, eye(n1)];
%     Ct = [Ct, ones(n1, n1 - n2)];
% end

% L
L = [KQ, -KQ * G2'; -G1 * KQ, G1 * KQ * G2' + KP];

% binary matrix -> index matrix (for saving memeory)
[IndG1, IndG1T, IndH1, IndH1T] = mat2inds(G1, G1', H1, H1');
[IndG2, IndG2T, IndH2, IndH2T] = mat2inds(G2, G2', H2, H2');
