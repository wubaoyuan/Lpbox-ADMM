function demo_normalizeMatchingW

% create a random graph matching matrix with n1 nodes in 1st graph, n2 in 2nd graph
n1 = 10;
n2 = 15;

% generate random feasible matches between 2 graphs
E12 = rand(n1, n2) > 0.3; % E12 = ones(n1, n2) for full matching
n12 = nnz(E12); % nb of feasible matches

% W(e, e') = affinity between match (i1,i2) and match (i1',i2')
% we only need to have a n12 x n12 matrix, as opposed to a 
% (n1 * n2) x (n1 * n2) matrix for memory efficiency and performance)
W = rand(n12, n12);
W = (W + W') / 2;

% assuming symmetric affinity, only need to keep upper triangular part
W = sparse(tril(W));
% in practice, memory-efficient code will directly generate a sparse
% triangular matrix

% kronecker normalization
normalization = 'iterative';
nbIter = 10;
[W, D1, D2] = normalizeMatchingW(W, E12, normalization, nbIter);

% verify that kronecker normalization worked (not really needed, just for the demo)
W = trilW2W(W); % inverse operation from tril(W)
[I12(:, 1), I12(:, 2)] = find(E12);
indexes1 = classes2indexes(I12(:, 1), n1); %indexes1{i1}=set of matches e=i1i2 in I12 (across i2)
indexes2 = classes2indexes(I12(:, 2), n2); %indexes2{i2}=set of matches e=i1i2 in I12 (across i1)
Dn1 = compute_D(W, indexes1); %Dn1(i1,j1)=sum of W(i1i2,j1j2) (across i2,j2) 
Dn2 = compute_D(W, indexes2); %Dn2(i2,j2)=sum of W(i1i2,j1j2) (across i1,j1) 

% verify that it is near-constant
std(Dn1(:))
std(Dn2(:))
0;
% [W,E12,target,G1,G2] = computeSyntheticExampleMatching(n1,n2,param,noise);

function D = compute_D(W, indexes)
n = length(indexes);
D = zeros(n);
for j = 1 : n
    for i = 1 : n
        D(i, j) = sum(sum(W(indexes{i}, indexes{j})));
    end
end
