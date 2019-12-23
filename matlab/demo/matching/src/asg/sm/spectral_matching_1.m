%% author: Marius Leordeanu
%  last updated: Feb 25, 2011

%  for questions contact the author at: leordeanu@gmail.com

%% please cite the following papers:
%  
%  A spectral technique for correspondence problems using pairwise
%  constraints (ICCV 2005)
%  by Marius Leordeanu ,  Martial Hebert 

% Utility:
% this function tries to maximize the matching score x'Mx 
% where x obeys discrete one-to-one matching constraints 
% such that x(i) = 1 if nodes(i) is matched to labels(i) and 0 otherwise

function [sol, v] = spectral_matching_1(M, labels, nodes)

minRatio = eps;
v = ones(length(nodes), 1);
v = v / norm(v);
iterClimb = 30;
nNodes = max(nodes);
nLabels = max(labels);

% compute the first eigenvector (iterative power method)
for i = 1 : iterClimb
    v = M * v;
    v = v / norm(v);
end

% make v double stochastic
aux = v;
v0 = aux;
v1 = aux;

for k = 1 : 20
    for j = 1 : nNodes
        f = find(nodes == j);
        v1(f) = v0(f) / (sum(v0(f)) + eps);
    end
    for j = 1 : nLabels
        f = find(labels == j);
        v0(f) = v1(f) / (sum(v1(f)) + eps);
    end
end
v = (v1 + v0) / 2;
v = v / norm(v);

% discretization using hungarian method 
A = zeros(nNodes, nLabels);
for i = 1 : nNodes
    f = find(nodes == i);
    A(i, :) = v(f);
end
A = max(A(:)) - A;
[X, score] = hungarian(A);
sol = X';
sol = sol(:);
