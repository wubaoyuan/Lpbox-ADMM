function vs = conRBF(X, Y, weis, sig, varargin)
% Compute (squared) distance matrix.
%
% Input
%   X       -  1st sample matrix, d x nX
%   Y       -  2nd sample matrix, d x nY
%   weis    -  weights, 1 x nY
%   sig     -  sigma
%   varargin
%     len   -  length of block, {4000}
%
% Output
%   vs      -  distance, nX x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-30-2012

% function option
len = ps(varargin, 'len', 4000);

% dimension
nX = size(X, 2);
nY = size(Y, 2);
XX = sum(X .* X);
YY = sum(Y .* Y); 

% divide
[idxXs, mX] = divCo2(nX, len);
[idxYs, mY] = divCo2(nY, len);

% compute each block
vs = zeros(nX, 1);
for iX = 1 : mX
    idxX = idxXs{iX};
    nXi = length(idxX);

    Xi = X(:, idxX);
    XXi = XX(:, idxX);
    for iY = 1 : mY
        idxY = idxYs{iY};
        nYi = length(idxY);

        Yi = Y(:, idxY);
        YYi = YY(:, idxY);

        weiY = weis(idxY);

        XY = Xi' * Yi;
        D = repmat(XXi', [1, nYi]) + repmat(YYi, [nXi, 1]) - 2 * XY;
        K = exp(-D .* sig) .* repmat(weiY, nXi, 1);

        vs(idxX) = vs(idxX) + sum(K, 2);
    end
end
