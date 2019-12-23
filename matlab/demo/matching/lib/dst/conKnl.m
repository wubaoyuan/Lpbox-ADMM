function [K, sigma] = conKnl(D, varargin)
% Construct the kernel matrix from squared distance matrix.
%
% Input
%   D       -  squared distance matrix, n x n
%              If D == [], then a distance matrix has been put as a global variable.
%   varargin 
%     knl   -  kernel type, 'g' | 'st' | 'bi' | 'man'
%              'g':  Gaussian kernel
%              'st': Self Tunning kernel
%                    see the paper, Self-Tuning Spectral Clustering, in NIPS 2004
%              'bi': Bipartite Maximum Weight Matching
%                    see the paper, Loopy Belief Propagation for Bipartite Maximum Weight b-Matching, in AISTATS 2007
%              'man': manually specified
%     nei   -  #nearest neighbour to compute the kernel bandwidth, {.1}
%              NaN: set bandwidth to 1
%              0: binary kernel
%              see function "bandG" for more details
%
% Output
%   K       -  kernel matrix, n x n
%              If nargout == 0, then put the kernel matrix as a global variable.
%   sigma   -  kernel bandwidth
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-04-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-03-2012

% global
global DG KG;
isDG = isempty(D);
isKG = nargout == 0;

% function option
knl = ps(varargin, 'knl', 'g');
nei = ps(varargin, 'nei', .1);

% Gaussian kernel
if strcmp(knl, 'g')
    sigma = bandG(D, nei);
    if isDG && isKG
        KG = exp(- DG ./ (2 * sigma ^ 2 + eps));
    elseif isDG && ~isKG
        K = exp(- DG ./ (2 * sigma ^ 2 + eps));
    elseif ~isDG && isKG
        KG = exp(- D ./ (2 * sigma ^ 2 + eps));
    else
        K = exp(- D ./ (2 * sigma ^ 2 + eps));
    end
    
% Gaussian kernel
elseif strcmp(knl, 'man')
    sigma = nei;
    if isDG && isKG
        KG = exp(- DG ./ (2 * sigma ^ 2 + eps));
    elseif isDG && ~isKG
        K = exp(- DG ./ (2 * sigma ^ 2 + eps));
    elseif ~isDG && isKG
        KG = exp(- D ./ (2 * sigma ^ 2 + eps));
    else
        K = exp(- D ./ (2 * sigma ^ 2 + eps));
    end

elseif strcmp(knl, 'st')
    Dsorted = sort(D);

    D2 = Dsorted(1 : m, :);
    sigma = sum(D2, 1) / m;

    K = exp(- D ./ (sigma' * sigma * 2 + eps));

elseif strcmp(knl, 'bi')
    lb = m;
    ub = m;
    mask = ones(n, n);
    P = bdmatchwrapper(-D, mask, lb, ub);

    D2 = D .* P;
    sigma = sum(D2, 1) / m;
    K = exp(- (D .^ 2) ./ (sigma' * sigma * 2 + eps));

else
    error('unknown kernel type: %s', knl);
end
