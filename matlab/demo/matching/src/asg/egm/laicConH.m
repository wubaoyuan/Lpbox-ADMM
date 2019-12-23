function [Hs, weis, Err] = laicConH(gphs, Use, varargin)
% Compute the reconstruciton matrix for LAIC algorithm.
%
% Input
%   gphs    -  1 x m (cell)
%   Use     -  flag of used graphs, k x m
%   varargin
%     eta   -  eta for reguralization, {1000}
%
% Output
%   Hs      -  reconstruction matrix, 1 x 1 (cell), k x k
%   weis    -  weights, 1 x k
%   Err     -  error, k x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-22-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-04-2012

% function option
eta = ps(varargin, 'eta', 1000);

% dimension
m = length(gphs);
k = size(gphs{1}.Pt, 2);

% default Use
if isempty(Use)
    Use = ones(k, m);
end

% template graph
[Pts, As] = cellss(1, m);
for i = 1 : m
    [Pt, Eg] = stFld(gphs{i}, 'Pt', 'Eg');

    % edge -> adjacency matrix
    A0 = gphEg2Adj(Eg, k);
    
    % make sure every node has at least three neighbours
    A = laicValA(A0);
    
    Pts{i} = Pt;
    As{i} = A;
end

% calculate the reconstruction matrix
[H, weis, Err] = calcReconCoes(Pts, As, Use, eta);
Hs = {H};
weis = weis';
