function Y = kpca(K, varargin)
% Kernel Principal Componant Analysis (KPCA).
%
% Input
%   K       -  kernel matrix, n x n
%   varargin
%     cen   -  flag of the fact that the input need to be centralized, {'y'} | 'n'
%     sel   -  method of selecting component, {'top'} | 'egy'
%     top   -  #top components
%     egy   -  energy to preserve
%
% Output
%   Y       -  new sample matrix after embedding, m x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-30-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
isCen = psY(varargin, 'cen', 'y');
sel = ps(varargin, 'sel', 'top');

% centralize kernel matrix
if isCen
    K = cenK(K);
end

% top component
if strcmp(sel, 'top')
    b = ps(varargin, 'top', 2);    
    V = eigk(K, b);

% preserve certain energy
elseif strcmp(sel, 'egy')
    [V, D] = eig(K);

    % sort
    [Lamb, index] = sort(sum(D, 2), 'descend');
    V = V(:, index);

    b = thEgy(Lamb, egy);

else
    error('unknown selecting method');
end

Y = V(:, 1 : b)' * K;
