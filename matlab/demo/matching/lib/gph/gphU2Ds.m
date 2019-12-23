function gphDs = gphU2Ds(gphUs)
% Convert an undirected graph to a directed one.
%
% Input
%   gphUs   -  undirected graphs, 1 x m (cell)
%
% Output
%   gphDs   -  directed graphs, 1 x m (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

% dimensin
m = length(gphUs);

gphDs = cell(1, m);
for i = 1 : m
    gphDs{i} = gphU2D(gphUs{i});
end
