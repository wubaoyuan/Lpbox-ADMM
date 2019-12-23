function [haPt, haEg] = plotGph(Pt, Eg, parMk)
% Plot one graph in 2-D space.
%
% Input
%   Pt       -  graph node, 2 x n
%   Eg       -  graph edge, 2 x m | []
%   parMk    -  marker parameter
%     cl     -  edge color, {'r'}
%     mk     -  marker type, {'o'} | 's' | ...
%     mkSiz  -  marker size, {5}
%     lnWid  -  line width, {1}
%            
% Output     
%   haPt     -  handle for nodes
%   haEg     -  handle for edge
%            
% History    
%   create   -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 01-02-2013

% function parameter
cl = ps(parMk, 'cl', 'r');
egCl = ps(parMk, 'egCl', 'r');
mkCl = ps(parMk, 'mkCl', 'r');
mk = ps(parMk, 'mk', 'o');
mkSiz = ps(parMk, 'mkSiz', 5);
lnWid = ps(parMk, 'lnWid', 1);

% dimension
m = size(Eg, 2);

hold on;

% plot edge
haEg = [];
if m > 0 && lnWid > 0
    X = [Pt(1, Eg(1, :)); Pt(1, Eg(2, :)); nan(1, m)];
    Y = [Pt(2, Eg(1, :)); Pt(2, Eg(2, :)); nan(1, m)];
    haEg = plot(X(:), Y(:), '-', 'color', cl, 'LineWidth', lnWid);
end

% plot node
haPt = [];
if mkSiz > 0
    haPt = plot(Pt(1, :), Pt(2, :), mk, 'MarkerFaceColor', mkCl, ...
                'MarkerEdgeColor', egCl, 'MarkerSize', mkSiz);
end
