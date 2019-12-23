function [Pts, Egs, asgs] = gphReOrd(Pt0s, Eg0s, asg0s)
% Re-order the nodes of each graph based on 1st graph's 2-D position.
% 
% Input
%   Pt0s    -  original points, 1 x 2 (cell)
%   Eg0s    -  original edges, 1 x 2 (cell)
%   asg0s   -  assignment, 1 x mG (cell) | struct
%
% Output
%   Pts     -  new points, 1 x 2 (cell)
%   Egs     -  new edges, 1 x 2 (cell)
%   asgs    -  assignment
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-19-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-20-2012

Pts = Pt0s;
Egs = Eg0s;
asgs = asg0s;

% distances
dists = sum(Pts{1} .^ 2, 1);
[~, ord] = sort(dists);

% adjust the nodes
Pts{1} = Pts{1}(:, ord);

% adjust the edges
if ~isempty(Egs{1})
    Pt = Pts{1};
    Eg = Egs{1};
    
    % edge -> adjacency
    n = size(Pt, 2);
    A = zeros(n, n);
    idx = sub2ind([n n], Eg(1, :), Eg(2, :));
    A(idx) = 1;
    A = triu(A, 1);
    A = A + A';
    
    % re-order
    A = A(ord, ord);
    
    % adjacency -> edge
    idx = find(A);
    [is, js] = ind2sub([n n], idx);
    Egs{1} = [is'; js'];
end

% adjust the correspondence
if iscell(asgs)
    for i = 1 : length(asgs)
        asgs{i}.X = asgs{i}.X(ord, :);
    end
else
    asgs.X = asgs.X(ord, :);    
end

