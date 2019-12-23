function ha = plotmk(X, c, parMk)
% Plot points belonging to one cluster in 2-D or 3-D.
%
% Input
%   X        -  sample matrix, d (= 2 | 3) x n
%   c        -  class label, 1 x 1
%   parMk    -  marker parameter
%     ln     -  line, {[]}
%               If empty, using lns{c} instead
%     mk     -  marker, {[]}
%               If empty, using mks{c} instead
%     cl     -  color, {[]}
%               If empty, using cls{c} instead
%     lns    -  line, {'-', '--', ':', '-.'}
%     mks    -  makrer, {'o', 's', '^', 'd', '+', 'v', 'x', '*', 'p'}
%     cls    -  color,  {[1 0 0], [0 0 1], [0 1 0], [1 0 1], [0 0 0], [0 1 1], [.3 .3 .3], [.5 .5 .5], [.7 .7 .7], [.1 .1 .1], [1 .8 0], [1, .4, .6]}
%     lnWid  -  line width, {1}
%     mkSiz  -  marker size, {5}
%     mkEg   -  marker edge, 'y' | {'n'}
%     face   -  flag of showing face, {'y'} | 'n'
%
% Output
%   ha       -  plot handle
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 11-19-2010
%   modify   -  Feng Zhou (zhfe99@gmail.com), 06-24-2012

% function parameter
ln = ps(parMk, 'ln', []);
mk = ps(parMk, 'mk', []);
cl = ps(parMk, 'cl', []);
lns = ps(parMk, 'lns', {'-', '--', ':', '-.'});
mks = ps(parMk, 'mks', {'o', 's', '^', 'd', '+', 'x', 'v', '>', '*'});
cls = ps(parMk, 'cls', {[1 0 0], [0 0 1], [0 1 0], [1 0 1], [0 0 0], [1 .5 0], [.7 .7 .7], [.1 .1 .1], [.4 .4 .7], [.1 .1 .1], [1 .8 0], [1, .4, .6]});
lnWid = ps(parMk, 'lnWid', 1);
mkSiz = ps(parMk, 'mkSiz', 5);
mkEgWid = ps(parMk, 'mkEgWid', 0);
mkEgCl = ps(parMk, 'mkEgCl', 'k');
isFace = psY(parMk, 'face', 'y');

% line
if isempty(ln)
    cc = mod(c - 1, length(lns)) + 1;
    ln = lns{cc};
end

% marker
if isempty(mk)
    cc = mod(c - 1, length(mks)) + 1;
    mk = mks{cc};
end

% color
if isempty(cl)
    cc = mod(c - 1, length(cls)) + 1;
    cl = cls{cc};
end

% dimension
d = size(X, 1);
if d ~= 2 && d ~= 3
    error('not supported dimension: d %d', d);
end

% plot line and marker
if lnWid > 0 && mkSiz > 0
    if d == 2
        ha = plot(X(1, :), X(2, :), ln);
    else
        ha = plot3(X(1, :), X(2, :), X(3, :), ln);
    end
    if isFace
        set(ha, 'Color', cl, 'LineWidth', lnWid, 'Marker', mk, 'MarkerSize', mkSiz, 'MarkerFaceColor', cl);
    else
        set(ha, 'Color', cl, 'LineWidth', lnWid, 'Marker', mk, 'MarkerSize', mkSiz);
    end
    if mkEgWid > 0        
        set(ha, 'MarkerEdgeColor', 'k');
    end

% plot marker only
elseif lnWid == 0 && mkSiz > 0
    if d == 2
        ha = plot(X(1, :), X(2, :), mk);
    else
        ha = plot3(X(1, :), X(2, :), X(3, :), mk);
    end
    
    % face
    if isFace
        set(ha, 'Color', cl, 'Marker', mk, 'MarkerSize', mkSiz, 'MarkerFaceColor', cl);
    else
        set(ha, 'Color', cl, 'Marker', mk, 'MarkerSize', mkSiz);
    end
    
    % marker edge
    if mkEgWid > 0
        set(ha, 'MarkerEdgeColor', mkEgCl, 'linewidth', mkEgWid);
    end

% plot line only
elseif lnWid > 0 && mkSiz == 0
    if d == 2
        ha = plot(X(1, :), X(2, :), ln);
    else
        ha = plot3(X(1, :), X(2, :), X(3, :), ln);
    end
    set(ha, 'Color', cl, 'LineWidth', lnWid);

else
    ha = [];
%    error('unsupported');
end
