function N = initNeighbors(I)
%%%%%%%%%%%%%%%%%%%%
% initNeighbors.m
% make a table of whose neighbors
%%%%%%%%%%%%%%%%%%%%
N = cell(numel(I), 1);
[r, c] = size(I);
for i = 1: numel(I)
    N{i} = getNeighbors(i, r, c);
end


