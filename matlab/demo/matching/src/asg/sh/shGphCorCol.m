function ha = shGphCorCol(gphs, X, XT, parCor)
% Show correspondence by different colors.
%
% Remark
%   n = n1 x n2
%
% Input
%   gphs     -  graph, 1 x 2 (cell)
%     Pt     -  graph node, 2 x ni
%   X        -  correspondence, n1 x n2
%   XT       -  ground-truth correspondence, [] | n1 x n2
%   parCor   -  parameter
%     mkSiz  -  marker size, {5}
%            
% Output     
%   ha       
%     HCor   -  handle of correspondence line, k10 x k20 (cell)
%     cls    -  color for line, 1 x 2 (cell)
%     X      -  correspondence, k10 x k20
%     XT     -  ground-truth correspondence, k10 x k20
%            
% History    
%   create   -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 02-15-2012

% function option
mkSiz = ps(parCor, 'mkSiz', 5);

% dimension
[n1, n2] = size(X);

% default ground-truth
if isempty(XT)
    XT = zeros(n1, n2);
end

% #correspondence
k = length(find((X + XT) > 0));
Cl = jet(k);

% per pair
HCor = cell(n1, n2);
co = 0;
for i = 1 : n1
    for j = 1 : n2
        % no correspondence
        if X(i, j) == 0 && XT(i, j) == 0
            continue;
        end
        co = co + 1;

        % coordinate
        pt1 = gphs{1}.Pt(:, i);
        pt2 = gphs{2}.Pt(:, j);
        cl = Cl(co, :);

        % plot
        plot(pt1(1), pt1(2), 'o', 'Color', cl, 'MarkerFaceColor', cl, 'MarkerSize', mkSiz, 'MarkerEdgeColor', cl);
        plot(pt2(1), pt2(2), 'o', 'Color', cl, 'MarkerFaceColor', cl, 'MarkerSize', mkSiz, 'MarkerEdgeColor', cl);
    end
end

% store
ha.HCor = HCor;
ha.X = X;
ha.XT = XT;
