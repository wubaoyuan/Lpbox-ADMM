function gphUs = gphD2Us(gphDs)
% Convert a directed graph to an undirected one.
%
% Input
%   gphDs   -  directed graphs, 1 x m (cell)
%
% Output
%   gphUs   -  undirected graphs, 1 x m (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% dimensin
m = length(gphDs);

gphUs = cell(1, m);
for i = 1 : m
    gphUs{i} = gphD2U(gphDs{i});
end
