function ha = shPtCol(Pt, lPt, k, parCor)
% Show point correspondence indicated with different colors.
%
% Input
%   Pt       -  graph node, 2 x n
%   lPt      -  node order, 1 x n
%   parCor   -  parameter
%     mkSiz  -  marker size, {5}
%            
% History    
%   create   -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 01-18-2012

% function option
mkSiz0 = ps(parCor, 'mkSiz', 5);

% dimension
n = size(Pt, 2);

% #correspondence
Cl = jet(k);

% per point
haPts = cell(1, n);
for i = 1 : n
    c = lPt(i);

    % coordinate
    pt = Pt(:, i);
    
    % color & marker
    if c == 0
        cl = [.5 .5 .5];
        mkSiz = round(mkSiz0 / 2);
    else
        cl = Cl(c, :);
        mkSiz = mkSiz0;
    end
    
    % plot
    haPts{i} = plot(pt(1), pt(2), 'o', 'Color', cl, 'MarkerFaceColor', cl, 'MarkerSize', mkSiz, 'MarkerEdgeColor', cl);
%    plot(pt(1), pt(2), 'o', 'Color', cl, 'MarkerFaceColor', cl, 'MarkerSize', mkSiz, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
end

% store
ha.haPts = haPts;
