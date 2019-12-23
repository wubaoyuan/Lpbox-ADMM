function [W, I12] = compute_matchingW(G1, W1, G2, W2, F12, W12, f, options)
% options.cliqueModeString={"log_linear","L1","edges","binary"}
% options.normalization={"none","D1","D2","D1D2","iterative"}
% options.nbIter=int
% f.fG=double[size(G,3)]
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

% TODO : arrange fact that G1,G2,F12 could be sparse

W12 = double(full(W12));
[I12(:, 1), I12(:, 2)] = find(W12);
n12 = sum(W12(:) > 0);
W12(W12(:) > 0) = 1 : n12;

% symmetric
n1 = length(W1);
W1 = sparse(W1 + W1' + speye(n1));
n2 = length(W2);
W2 = sparse(W2 + W2' + speye(n2));

G1 = full(double(G1));
G2 = full(double(G2));
F12 = full(F12);

if nargin < 7 || length(f.fG) ~= size(G1, 3)
    error('bad f');
end

W = mex_matchingW2(G1, W1, G2, W2, F12, I12, W12, f, options);

% [indi,indj] = find(W);
% ind = find(W);
% noise = 1e-6;
% noise = 1e-9;
% % W(ind) = W(ind) + noise*rand(length(ind),1);% to stabilize ???
% W = W + sparse(indi,indj,noise*rand(length(indi),1),n12,n12);% to stabilize ???


% v1 = ones(n1,1);
% v2 = ones(n2,1);
% r1=eye(n1)-1/n1*v1*v1';
% r2=eye(n2)-1/n2*v2*v2';
% R=kron(r2,r1);
% W = trilW2W(W);
% W=sparse(R*full(W)*R);
