function shHstCur(Hst, varargin)
% Plot histogram as curves in 2-D.
%
% Input
%   Hst       -  histogram, nRep x n x nAlg
%   varargin
%     show option
%     algs    -  algorithm names, {[]}
%     xs      -  x axis value, {[]}
%     cs      -  class of algorithm, {[]}
%     lnWid   -  line width, {1}
%     lns     -  line markers, {[]}
%     mkSiz   -  size of markers, {5}
%     mkEg    -  flag of marker edge, 'y' | {'n'}
%     Var     -  variances, {[]}
%     varWid  -  variance width in terms of axis gap, {.1}
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 01-15-2010
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% show option
psSh(varargin);

% function option
algs = ps(varargin, 'algs', []);
xs = ps(varargin, 'xs', []);
cs = ps(varargin, 'cs', []);
lnWid = ps(varargin, 'lnWid', 1);
mkSiz = ps(varargin, 'mkSiz', 5);
isMkEg = psY(varargin, 'mkEg', 'n');
varWid = ps(varargin, 'varWid', .1);
plFun = ps(varargin, 'plFun', 'plot');
plFun2 = ps(varargin, 'plFun2', 'plot');

% dimension
[nRep, n, nAlg] = size(Hst);

% xs
if isempty(xs)
    xs = 1 : n;
end

% cs
if isempty(cs)
    cs = 1 : nAlg;
end

% plot mean
mes = cell(1, nAlg);
for iAlg = 1 : nAlg
    Hsti = Hst(:, :, iAlg);
    vis = isnan(Hsti(1, :));
    Hsti(:, vis) = [];

    xi = xs(~vis);
    me = mean(Hsti, 1);
    mes{iAlg} = me;
    c = cs(iAlg);

    if strcmp(plFun, 'plot')
        hTmp = loglog(xi, me, '-', 'LineWidth', lnWid);
    elseif strcmp(plFun, 'loglog')
        hTmp = loglog(xi, me, '-', 'LineWidth', lnWid);
    elseif strcmp(plFun, 'semilogx')
        hTmp = semilogx(xi, me, '-', 'LineWidth', lnWid);
    elseif strcmp(plFun, 'semilogy')
        hTmp = semilogy(xi, me, '-', 'LineWidth', lnWid);        
    else
        error('unknown plot fun: %s', plFun);
    end

    % color
    [mk, cl] = genMkCl(c);
    set(hTmp, 'Color', cl);
    if mkSiz > 0
        set(hTmp, 'Marker', mk, 'MarkerSize', mkSiz, 'MarkerFaceColor', cl);
        if isMkEg
            set(hTmp, 'MarkerEdgeColor', 'k');
        end
    end
    
    if iAlg == 1
        hold on;
    end
end

% legend
if ~isempty(algs)
    legend(algs{:});
end
%return;

% plot variance
dx = diff(xs);
wid2 = dx(1);
for iAlg = 1 : nAlg
    Hsti = Hst(:, :, iAlg);
    vis = isnan(Hsti(1, :));
    Hsti(:, vis) = [];

    xi = xs(~vis);
    de = std(Hsti, 0, 1);
    me = mes{iAlg};
    c = cs(iAlg);
    ni = length(xi);

    % color
    [~, cl] = genMkCl(c);
    for i = 1 : ni
        if strcmp(plFun2, 'plot')
            plot([xi(i) xi(i)], [-de(i) de(i)] + me(i), '-', 'Color', cl, 'LineWidth', lnWid);
            plot([-wid2, wid2] * varWid + xi(i),  [1 1] * de(i) + me(i), '-', 'Color', cl, 'LineWidth', lnWid);
            plot([-wid2, wid2] * varWid + xi(i), -[1 1] * de(i) + me(i), '-', 'Color', cl, 'LineWidth', lnWid);
            
        elseif strcmp(plFun2, 'loglog')
            loglog([xi(i) xi(i)], [-de(i) de(i)] + me(i), '-', 'Color', cl, 'LineWidth', lnWid);
            loglog([-wid2, wid2] * varWid + xi(i),  [1 1] * de(i) + me(i), '-', 'Color', cl, 'LineWidth', lnWid);
            loglog([-wid2, wid2] * varWid + xi(i), -[1 1] * de(i) + me(i), '-', 'Color', cl, 'LineWidth', lnWid);
            
        elseif strcmp(plFun2, 'semilogx')
            semilogx([xi(i) xi(i)], [-de(i) de(i)] + me(i), '-', 'Color', cl, 'LineWidth', lnWid);
            semilogx([-wid2, wid2] * varWid + xi(i),  [1 1] * de(i) + me(i), '-', 'Color', cl, 'LineWidth', lnWid);
            semilogx([-wid2, wid2] * varWid + xi(i), -[1 1] * de(i) + me(i), '-', 'Color', cl, 'LineWidth', lnWid);
        end
    end
end
