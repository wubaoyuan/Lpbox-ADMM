function conKnlLi(X, dims, paraFt)
% Construct multiple kernel matrices from different dimensions of sample.
%
% Input
%   X       -  sample matrix, dim x n
%   dims    -  feature dimensions, 1 x nFt
%   paraFt  -  feature parameter
%     nei   -  #nearest neighbour to compute the kernel bandwidth, {.1}
%     neis  -  'nei' value for each feature
%     a     -  weight vector, 1 x nFt
%     Ks    -  flag of storing Ks, 'y' | {'n'}
%     nor   -  normalization flag, 'y' | {'n'}
%
% Output
%   K       -  kernel matrix, n x n
%   Ks      -  kernel matrices, n x n x nFt
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% global
global KG KGs;

% function option
neis = ps(paraFt, 'neis', []);
nei = ps(paraFt, 'nei', .1);
a = ps(paraFt, 'a', []);
isKs = psY(paraFt, 'Ks', 'n');
isNor = psY(paraFt, 'nor', 'n');

% data
[dim, n] = size(X);

% feature dimensions
if isempty(dims)
    dims = ones(1, dim);
end
nFt = length(dims);
s = n2s(dims);

% neighbor
if isempty(neis)
    neis = repmat(nei, nFt, 1);
end

% kernel of each part
if isKs
    KGs = zeross(n, n, nFt);
end

% X2s = cell(1, nFt);
KG = zeros(n, n);
for iFt = 1 : nFt
    idx = s(iFt) : s(iFt + 1) - 1;

    % distance
    conDst(X(idx, :), X(idx, :));

    % kernel
    Ki = conKnl([], 'nei', neis(iFt));

    if isKs
        KGs(:, :, iFt) = Ki;
    end
    KG = KG + Ki * a(iFt);

    % X
%     nn = sqrt(sum(X(idx, :) .^ 2));
%     X2s{iFt} = X(idx, :) ./ repmat(nn, length(idx), 1);
end
% X2 = cat(1, X2s{:});

% normalization
if isNor
    de2 = sqrt(sum(KG, 2));
    KG = KG ./ (de2 * de2');
end
