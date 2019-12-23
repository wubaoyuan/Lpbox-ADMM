function X = mapHist(X0, hist, varargin)
% Generate value based on the histogram.
%
% Input
%   X0      -  original sample matrix, dim x n
%   hist    -  histogram
%   varargin
%
% Output
%   X       -  new sample matrix, nBin (dim x val) x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-11-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
alg = ps(hist, 'alg', 'eq');
val = ps(hist, 'val', []);

[dim, n] = size(X0);

if strcmp(alg, 'eq')
    mi = hist.mi;
    ma = hist.ma;

    nBin = val * dim;
    L = zeros(dim, n);

    for di = 1 : dim
        maD = ma(di);
        miD = mi(di);
        
        if maD - miD < eps
            L(di, :) = ones(1, n);

        else
            for i = 1 : n
                L(di, i) = float2block((X0(di, i) - miD) / (maD - miD), val);
            end
        end
    end
    
elseif strcmp(alg, 'kmeans')
else
    error('unknown method');
end

X = zeros(nBin, 1);
for di = 1 : dim
    for i = 1 : val
        X((di - 1) * val + i) = length(find(L(di, :) == i));
    end
end
