function ha = shGphUpd(ha, gphs, varargin)
% Update multiple graphs in 2-D space.
%
% Input
%   gphs      -  graphs, 1 x k (cell)
%     Pt      -  graph node, 2 x ni
%     Eg      -  graph edge, 2 x mi | []
%   varargin
%     show option
%     parMks  -  marker parameter, {[]} | 1 x k (cell)
%     parAx   -  axis parameter, {[]}, see function setAx for more details
%     num     -  flag of whether ploting number, 'y' | {'n'}
%
% Output
%   ha        -  figure handle
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 02-15-2012

% function option
parMks = ha.parMks;
lPts = ps(varargin, 'lPts', []);

% dimension
k = length(gphs);

% default edge
if isempty(Egs)
    Egs = cell(1, k);
end

% update each graph
for i = 1 : k
    % update point's position
    set(haPts{i}, 'XData', Pts{i}(1, :), 'YData', Pts{i}(2, :));
end

% store
ha.haPts = haPts;
ha.haEgs = haEgs;
