function [G, cost, acc] = kmean(X, k, varargin)
% K-means. A wrapper of matlab function, kmeans.
% After initalized serveral times, the one with the minimum cost is selected.
%
% Input
%   X       -  sample matrix, dim x n
%   k       -  cluster number
%   varargin
%     nRep  -  number of repetition to select the minimum cost, {10}
%     GT    -  ground truth label, {[]} |  k x n
%
% Output
%   G       -  class indicator matrix, k x n
%   cost    -  error cost
%   acc     -  accuracy if GT is indicated
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-22-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
nRep = ps(varargin, 'nRep', 10);
GT = ps(varargin, 'GT', []);

XTran = X';
Gs = cell(1, nRep); 
costs = zeros(1, nRep);
warning off;
for i = 1 : nRep
    try
        L = kmeans(XTran, k, 'emptyaction', 'singleton', 'display', 'off');
    catch
        err = lasterror;
        sprintf('%s\n', err.message);
        
        costs(i) = inf;
        continue;
    end

    G = L2G(L, k);
    
    % cost
    costs(i) = evalClu(X, G, 'type', 'in');
    
    % store
    Gs{i} = G;
end
warning on;

% pick the one with least cost
[cost, ind] = min(costs);
G = Gs{ind};
acc = 0;

% ground-truth
if ~isempty(GT)
    [G, acc] = matchG(G, GT);
end
