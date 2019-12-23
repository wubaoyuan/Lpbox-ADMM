function [X, scores] = compute_ICM_graph_matching(X0, E12, W);
% maximum enrgy
% TODO:verifier formule; esp case when Ws_{ii}(xi,xi) full matrix???

% assert(~mex_istril(W));%TODO

if issparse(W)
    if(mex_istril(W))
        W = trilW2W(W);
        % W=(W+W')/2;
    end
end

[n, k] = size(X0);
if n > k
    [X, scores] = compute_ICM_graph_matching_transpose(X0, E12, W);
else
    [X, scores] = compute_ICM_graph_matching_direct(X0, E12, W);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, scores] = compute_ICM_graph_matching_transpose(X0, E12, W)
[n1, n2] = size(E12);
[I12(:, 1), I12(:, 2)] = find(E12);
[ignore, perm] = sort(I12(:, 2) + n2 * (I12(:, 1) - 1));
X0 = X0';
E12 = E12';

W = W(perm, perm);

[X, scores] = compute_ICM_graph_matching_direct(X0, E12, W);
X = X';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, scores] = compute_ICM_graph_matching_direct(X0, E12, W)

[n, k] = size(X0);
assert(n <= k);%TODO:support otherwise
scores = [];
X = X0(:);
D = diag(W);
D = full(D);
if issparse(W)
    Ws = W - spdiag(D);
else
    Ws = W - diag(D);
end

X = reshape(X, n, k);
assert(all(sum(X > 0, 2) == 1));
X0 = X;
[X, scores] = mex_ICM_graph_matching(X, E12, Ws, D);
