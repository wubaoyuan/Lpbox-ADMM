function X = norX2(X0, alg)
% Normalize each dimension spearately with respect to the variance of distance.
%
% Input
%   X0      -  original sample matrix, d x n
%   alg     -  algorithm, 'l1' | 'l2'
%
% Output
%   X       -  new sample matrix, d x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-22-2012

% dimension
[d, n] = size(X0);
len = 128;

if strcmp(alg, 'l1')
    % l1 normalization
    vs = sum(X0, 1) + eps;
    X = X0 ./ repmat(vs, d, 1);

elseif strcmp(alg, 'l12')
    X = zeros(d, n);
    m = round(d / len);
    
    % per bin
    for i = 1 : m
        pHd = (i - 1) * len + 1;
        pEd = i * len;
        idx = pHd : pEd;
        Xi = X0(pHd : pEd, :);
        
        % l2 normalization
        vs = sum(Xi .^ 2, 1) + eps;
        vs = real(sqrt(vs));
        Xi = Xi ./ repmat(vs, len, 1);
        
        % l1 normalization
        vs = sum(Xi, 1) + eps;
        Xi = Xi ./ repmat(vs, len, 1);
        
        X(idx, :) = Xi;
    end    

elseif strcmp(alg, 'l2')
    
else
    error('unknown alg, %s', alg);
end
