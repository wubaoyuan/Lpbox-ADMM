function boxs = gphBox(gphs, varargin)
% Obtain bounding box of graph.
% 
% Input
%   gphs    -  graph, 1 x m (cell)
%     Pt    -  graph nodes, 2 x ki
%   varargin
%
% Output
%   boxs    -  bounding box, 1 x m (cell), 2 x 2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-02-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-15-2012

% dimension
m = length(gphs);

% per graph
boxs = cell(1, m);
for i = 1 : m
    Pti = gphs{i}.Pt;
    boxs{i} = [min(Pti, [], 2), max(Pti, [], 2)];
end
