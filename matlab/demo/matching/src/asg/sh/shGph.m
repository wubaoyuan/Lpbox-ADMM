function ha = shGph(gphs, varargin)
% Show multiple graphs in 2-D space.
%
% Input
%   gphs      -  graphs, 1 x k (cell)
%     Pt      -  graph node, 2 x ni
%     Eg      -  graph edge, 2 x mi | []
%   varargin
%     show option
%     parMks  -  marker parameter, {[]} | 1 x k (cell)
%     parAx   -  axis parameter, {[]}, see function setAx for more details
%
% Output
%   ha        -  figure handle
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 02-16-2012

% show option
psSh(varargin);

% function option
parMks = ps(varargin, 'parMks', []);
parAx = ps(varargin, 'parAx', []);

% dimension
k = length(gphs);

% node & edge
[Pts, Egs] = cellss(1, k);
for c = 1 : k
    Pts{c} = gphs{c}.Pt;
    Egs{c} = gphs{c}.Eg;
end

% default marker parameter
if isempty(parMks)
    parMks = cell(1, k);
    for i = 1 : k
        parMks{i} = st('mkSiz', 5);
    end
end

% default axis parameter
if isempty(parAx)
    parAx = st('set', 'n');
end

hold on;

% plot each graph
[haPts, haEgs] = cellss(1, k);
for i = 1 : k
    [haPts{i}, haEgs{i}] = plotGph(Pts{i}, Egs{i}, parMks{i});
end

% axis
Pt = mcat('horz', Pts);
box = xBox(Pt, parAx);
setAx(box, parAx);

% store
ha.haPts = haPts;
ha.haEgs = haEgs;
