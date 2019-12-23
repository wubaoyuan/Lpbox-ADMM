function gphs = gphBoxNew(gph0s, box0s, boxs, Sca)
% Obtain new graph.
% 
% Input
%   gph0s   -  original graph, 1 x m (cell)
%     Pt    -  graph node, 2 x ki
%   box0s   -  original box set, 1 x m (cell), 2 x 2
%   boxs    -  new box set, 1 x m (cell), 2 x 2
%   Sca     -  scale (= siz / siz0), 2 x m
%
% Output
%   gphs    -  new graph, 1 x m (cell)
%     Pt    -  graph node, 2 x ki
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-02-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-15-2012

% dimension
m = length(gph0s);

% new node
gphs = gph0s;
for i = 1 : m
    gphs{i}.Pt = xBoxNew(gph0s{i}.Pt, box0s{i}, boxs{i}, Sca(:, i));
end
