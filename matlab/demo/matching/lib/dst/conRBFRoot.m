function vs = conRBFRoot(X, Y, weis, sig, varargin)
% Compute (squared) distance matrix for root SVM.
%
% Input
%   X       -  1st sample matrix, d x nX
%   Y       -  2nd sample matrix, d x nY
%   weis    -  weights, 1 x nY
%   sig     -  sigma
%   varargin
%     len   -  size of block, {4000}
%
% Output
%   vs      -  distance, nX x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-05-2012

% function option
len = ps(varargin, 'len', 4000);

% dimension
nX = size(X, 2);
nY = size(Y, 2);

% dividefa
[idxXs, mX] = divCo2(nX, len);
[idxYs, mY] = divCo2(nY, len);

% compute each block
vs = zeros(nX, 1);
for iX = 1 : mX
    idxX = idxXs{iX};
    nXi = length(idxX);

    Xi = X(:, idxX);
    for iY = 1 : mY
        idxY = idxYs{iY};
        Yi = Y(:, idxY);
        weiY = weis(idxY);

        XY = Xi' * Yi;
        D = 2 - 2 * XY;
        K = exp(-D .* sig) .* repmat(weiY, nXi, 1);

        vs(idxX) = vs(idxX) + sum(K, 2);
    end
end
