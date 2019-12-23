function [Hst, Idx] = ctx2hst(Ctxs, par)
% Compute histogram from shape context.
%
% Input
%   Ctxs     -  shape context, 1 x nF (cell), nBin x nPt
%   par      -  function parameter
%     k      -  #class for k-means, {20}
%     nItMa  -  #maximum iterations, {200}
%     nRep   -  #repetitions, {1}
%
% Output
%   Hst      -  histogram of shape context, k x nF
%   Idx      -  label index of each point, nPt x nF
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 09-26-2010
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function parameter
k = ps(par, 'k', 20);
nItMa = ps(par, 'nItMa', 200);
nRep = ps(par, 'nRep', 1);

% dimension
nF = length(Ctxs);
nPt = size(Ctxs{1}, 2);
Ctx = mcat('horz', Ctxs);

% build dictionary by k-means
opts = statset('maxiter', nItMa);
idx = kmeans(Ctx', k, 'emptyaction', 'singleton', 'display', 'off', ...
             'Replicates', nRep, 'options', opts);

% histogram
Idx = zeros(nPt, nF);
Hst = zeros(k, nF);
for iF = 1 : nF
    for iPt = 1 : nPt
        c = idx(nPt * (iF - 1) + iPt);
        Idx(iPt, iF) = c;
        Hst(c, iF) = Hst(c, iF) + 1;
    end
end
