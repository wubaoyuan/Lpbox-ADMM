function W = trilW2W(W)

if issparse(W)
    temp = spdiag(diag(W));
else
    temp = diag(diag(W));
end
W = W + W';
W = W - temp;
