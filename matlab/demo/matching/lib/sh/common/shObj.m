function shObj(objss, varargin)
% Show cost at each iterations.
%
% Input
%   objss    -  objective, 1 x m (cell), 1 x nIti
%   varargin
%     show option
%     algs   -  algorithm name, {[]}
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 09-20-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% show option
psSh(varargin);

% function option
algs = ps(varargin, 'algs', []);

% dimension
nIts = cellDim(objss, 2);

vis = nIts == 0;
objss(vis) = [];
nIts(vis) = [];
if ~isempty(algs)
    algs(vis) = [];
end

m = length(objss);

% plot
hold on;
for i = 1 : m    
    [mk, cl] = genMkCl(i);
    plot(1 : nIts(i), objss{i}, '-', 'Color', cl, 'Marker', mk);
end

% legend
if ~isempty(algs)
    legend(algs{:});
end
