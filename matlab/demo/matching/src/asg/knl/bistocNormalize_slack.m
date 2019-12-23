function [X, Xslack] = bistocNormalize_slack(X, tolC)
% Sinkhorn method of iterative bistocastic normalizations.
%

% dimension
[n1, n2] = size(X);

% non-square
if n1 ~= n2
    Xslack = X;
    if n1 > n2
        Xslack(:, n2 + 1 : n1) = 1;
    else
        Xslack(n1 + 1 : n2, :) = 1;
    end
    [Xslack, ~] = mex_normalize_bistochastic(Xslack, tolC, 1000);
    X = Xslack(1 : n1, 1 : n2);

% square    
else
    Xslack = X;
    
    [X, ~] = mex_normalize_bistochastic(X, tolC, 1000);
end

