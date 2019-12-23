function h = shAcc2(Acc, varargin)
% Plot histogram in 2-D space.
%
% Input
%   Accs      -  accuracy, nAlg x n
%   varargin
%     show option
%     leg     -  legend, {[]}
%     xs      -  x axis value, {[]}
%     cs      -  class of algorithm, {[]}
%     lnWid   -  line width, {1}
%     lns     -  line markers, {[]}
%     mkSiz   -  size of markers, {5}
%     mkEg    -  flag of marker edge, {'y'} | 'n'
%     Var     -  variances, {[]}
%     varWid  -  variance width in terms of axis gap, {.1}
%
% Output
%   h         -  handle
%     lns     -  lines, 1 x nAlg (cell)
%     leg     -  legend
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 01-15-2010
%   modify    -  Feng Zhou (zhfe99@gmail.com), 02-29-2012

% show option
psSh(varargin);

% function option
leg = ps(varargin, 'leg', []);
xs = ps(varargin, 'xs', []);
cs = ps(varargin, 'cs', []);
lnWid = ps(varargin, 'lnWid', 1);
mkSiz = ps(varargin, 'mkSiz', 5);
isMkEg = psY(varargin, 'mkEg', 'y');
Var = ps(varargin, 'Var', []);
varWid = ps(varargin, 'varWid', .1);
nG = ps(varargin, 'nG', []);

% dimension
[nAlg, n] = size(Acc);

% xs
if isempty(xs)
    xs = 1 : n;
end

% cs
if isempty(cs)
    cs = 1 : nAlg;
end

% markers
[markers, colors] = genMkCl;

% plot
hold on;
for iAlg = 1 : nAlg
    if isempty(nG)
        hTmp = plot(xs, Acc(iAlg, :), '-', 'LineWidth', lnWid);
    else
        if iAlg <= nG
            hTmp = plot(xs, Acc(iAlg, :), '-', 'LineWidth', lnWid);
        else
            hTmp = plot(xs, Acc(iAlg, :), '--', 'LineWidth', lnWid);
        end
    end
    
    % color
    if isempty(nG)
        c = cs(iAlg);
    else
        c = mod(iAlg - 1, nG) + 1;
    end
    cl = colors{1 + mod(c - 1, length(colors))};
    mk = markers{1 + mod(c - 1, length(markers))};

    if isempty(c)
        set(hTmp, 'Color', colors{1});
    else
        set(hTmp, 'Color', cl);

        % marker
        if mkSiz > 0
            set(hTmp, 'Marker', mk, 'MarkerSize', mkSiz, 'MarkerFaceColor', cl);

            if isMkEg
                set(hTmp, 'MarkerEdgeColor', 'k');
            end
        end
    end
    lns{iAlg} = hTmp;
end
h.lns = lns;

% legend
if ~isempty(leg)
    h.leg = legend(leg{:});
end

% variance
if ~isempty(Var)
    for iAlg = 1 : nAlg
        c = cs(iAlg);
        acc = Acc(iAlg, :);
        dev = sqrt(Var(iAlg, :));
        dx = diff(xs);
        wid2 = dx(1);
        for i = 1 : n
            line([xs(i) xs(i)], [-dev(i) dev(i)] + acc(i), 'Color', colors{c}, 'LineWidth', 1);
            line([-wid2, wid2] * varWid + xs(i), [1 1] * dev(i) + acc(i), 'Color', colors{c}, 'LineWidth', 1);
            line([-wid2, wid2] * varWid + xs(i), -[1 1] * dev(i) + acc(i), 'Color', colors{c}, 'LineWidth', 1);
        end
    end
end
