function h = shs(Xs, parMk, parAx, varargin)
% Show point set.
%
% Input
%   Xs      -  sequence set, 1 x m (cell), dim x ni
%   parMk   -  marker parameter
%   parAx   -  axis parameter
%   varargin
%     show option
%     G     -  class label, k x m
%
% Output
%   h       -  figure handle
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-31-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-04-2013

% show option
psSh(varargin);

% function option
G = ps(varargin, 'G', []);

% dimension
m = length(Xs);

% label
if isempty(G)
    G = eye(m);
end
l = G2L(G);

% main plot
hold on;
for i = 1 : m
    plotmk(Xs{i}, l(i), parMk);
end

% axis
h.box = xBox(cat(2, Xs{:}), parAx);
setAx(h.box, parAx);
