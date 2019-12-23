function [A, vis] = gphAdjDelRef(Pt0, nConn)
% Delaunay triangulation with refinement.
%
% Input
%   Pt0     -  original graph node, 2 x n0
%   nConn   -  minimum #connection number for each node
%
% Output
%   A       -  node-node adjacency, n x n
%   vis     -  binary indicator of nodes that have been kept, 1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-29-2012

prIn('gphAdjDel', 'nConn %d', nConn);
isDeb = 0;

% dimension
Pt = Pt0;
n0 = size(Pt, 2);
vis = ones(1, n0) > 0;

% remove points that are closer than threshold
D = conDst(Pt, Pt);
idx = sub2ind([n0 n0], 1 : n0, 1 : n0);
D(idx) = inf;
idx = find(D(:) < eps);
if length(idx) > 0
    error('some points are too closer');
end

% delete one node until satisfying the constraints
while true
    % current point
    Pt = Pt0(:, vis);
    n = size(Pt, 2);

    % core triangluation
    A = gphAdjDel(Pt);

    % #neighbors
    ks = sum(A, 2);

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
    if all(ks >= nConn)
        break;
    end

    % remove one point
    ii = find(ks < nConn, 1);
    idx = find(vis);
    vis(idx(ii)) = false;
end

prOut;