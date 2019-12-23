function ha = shCur(Me, Dev, varargin)
% Plot curve group in 2-D.
%
% Input
%   Me         -  mean matrix, k x n
%                   k : #curves (algorithms)
%                   n : #points (bins) per curve
%   Dev        -  standard variation matrix, k x n
%   varargin
%     show option
%     barWidG  -  width of group of histgram bar, {.8} <= 1
%     barWid   -  width of histgram bar in terms of 'barWidG / k', {.8} <= 1
%     devWid   -  width of deviation line in terms of the 'barWid', {1} <= 1
%     bdWid    -  width of boundary of bar, {1}
%     dev      -  dev flag, {'y'} | 'n'
%     yLim     -  threshold in y axis, {[]}
%     leg      -  legend, {[]}
%
% Output
%   ha         -  handle
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 03-03-2009
%   modify     -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

% show option
psSh(varargin);

% function option
devWid = ps(varargin, 'devWid', .1);
bdWid = ps(varargin, 'bdWid', 1);
parMk = ps(varargin, 'parMk', []);
isDev = psY(varargin, 'dev', 'n');
yLim = ps(varargin, 'yLim', []);
algs = ps(varargin, 'algs', []);
xs = ps(varargin, 'xs', []);

% body
hold on;

% default parMk
if isempty(parMk)
    parMk = st('ln', '-', 'mkSiz', 6, 'lnWid', 1);
end

% mk & cl
[mks, cls] = genMkCl;

% dimension
[k, n] = size(Me);

% width
devWid = devWid * 1;

% x position
if isempty(xs)
    xs = 1 : n;
end

may = max(max(Me));
miy = min(min(Me));

% standard deviation
if isDev
    cls = ps(parMk, 'cls', {[1 0 0], [0 0 1], [0 1 0], [1 0 1], [0 0 0], [1 .5 0], [.7 .7 .7], [.1 .1 .1], [.4 .4 .7], [.1 .1 .1], [1 .8 0], [1, .4, .6]});
    dx = diff(xs);
    wid2 = dx(1) / 2;
    for c = 1 : k
        for i = 1 : n
            plot([xs(i) xs(i)], [-1, 1] * Dev(c, i) + Me(c, i), 'Color', cls{c}, 'LineWidth', 1);
            plot([-wid2, wid2] * devWid + xs(i), [1 1] * Dev(c, i) + Me(c, i), 'Color', cls{c}, 'LineWidth', 1);
            plot([-wid2, wid2] * devWid + xs(i), -[1 1] * Dev(c, i) + Me(c, i), 'Color', cls{c}, 'LineWidth', 1);
        end
    end
    may = may + max(max(Dev));
end

% mean
for c = 1 : k
    plotmk([xs; Me(c, :)], c, parMk);
end

% legend
if ~isempty(algs)
    ha.leg = legend(algs{:}, 'orientation', 'horizontal');
end

% boundary
xlim([xs(1) - 1, xs(end) + 1]);
if isempty(yLim)
    gap = (may - miy) * .1;
    mi = miy - gap;
    ma = may + gap;
else
    mi = yLim(1); ma = yLim(2);
end
if abs(mi - ma) > eps
    ylim([mi, ma]);
end
