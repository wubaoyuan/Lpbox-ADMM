function [Pt, XP, Eg, A] = gphConn2Old(Pt0, XP0)
% Connect nodes to generate edges.
%
% Input
%   Pt0     -  original graph node, 2 x n0
%
% Output
%   Eg      -  graph edge, 2 x m
%   A       -  node-node adjacency, n x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-29-2012

prIn('gphConn2');
isDeb = 0;

% dimension
Pt = Pt0;
XP = XP0;
n = size(Pt, 2);

while true
    D = conDst(Pt, Pt);
    idx = sub2ind([n n], 1 : n, 1 : n);
    D(idx) = inf;
    idx = find(D(:) < eps);
    [is, js] = ind2sub([n n], idx);
    pr('length(idx): %d', length(idx));
    
    Tri = delaunay(Pt(1, :)', Pt(2, :)');
    nTri = size(Tri, 1);

    % adjacency matrix
    A = zeros(n, n);
    for iTri = 1 : nTri
        i1 = Tri(iTri, 1);
        i2 = Tri(iTri, 2);
        i3 = Tri(iTri, 3);
        is = sort([i1 i2 i3]);

        A(is(1), is(2)) = 1;
        A(is(1), is(3)) = 1;
        A(is(2), is(3)) = 1;
    end
    A = A + A';    
    
    % neighbors
    idxs = cellss(n, 1);
    for i = 1 : n
        idxs{i} = find(A(i, :) ~= 0);
    end
    ks = cellDim(idxs, 2);
    
    % plot
    if isDeb
        figure(101);
        clf;
        hold on;
        
        idx = find(A);
        [is, js] = ind2sub([n n], idx);
        Eg = [is, js]';
        for i = 1 : size(Eg, 2)
            plot(Pt(1, Eg(:, i)), Pt(2, Eg(:, i)), 'r-');
        end
        
        plot(Pt(1, ks >= 3), Pt(2, ks >= 3), 'ro');
        plot(Pt(1, ks < 3), Pt(2, ks < 3), 'bs', 'MarkerFaceColor', 'b');
        drawnow;
        pause(.1);
    end
    
    % stop condition
    if all(ks >= 3)
        break;
    end
    
    % remove point
    ii = find(ks < 3, 1);
    Pt(:, ii) = [];
    XP(:, ii) = [];
    n = size(Pt, 2);
end
A = triu(A, 1);

% index of edge
idx = find(A);
[is, js] = ind2sub([n n], idx);
Eg = [is, js]';

prOut;