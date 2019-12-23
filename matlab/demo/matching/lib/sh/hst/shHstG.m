function ha = shHstG(Me, Dev, varargin)
% Plot histogram group in 2-D.
%
% Input
%   Me         -  mean, k x n
%                   k : #bin in each group
%                   n : #group
%   Dev        -  deviation, k x n
%   varargin
%     show option
%     barWidG  -  width of group of histgram bar, {.8} <= 1
%     barWid   -  width of histgram bar in terms of 'barWidG / k', {.8} <= 1
%     devWid   -  width of deviation line in terms of the 'barWid', {1} <= 1
%     bdWid    -  width of boundary of bar, {1}
%     lnWid    -  width of deviation line, {2}
%     dev      -  dev flag, {'y'} | 'n'
%     val      -  value flag, 'y' | {'n'}
%     form     -  value form, {'%d'}
%     yLim     -  threshold in y axis, {[]}
%     leg      -  legend, {[]}
%
% Output
%   ha         -  handle
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 03-03-2009
%   modify     -  Feng Zhou (zhfe99@gmail.com), 04-06-2013

% show option
psSh(varargin);

% function option
barWidG = ps(varargin, 'barWidG', .8);
barWid = ps(varargin, 'barWid', .8);    
devWid = ps(varargin, 'devWid', .6);
bdWid = ps(varargin, 'bdWid', 1);
lnWid = ps(varargin, 'lnWid', 2);
isDev = psY(varargin, 'dev', 'n');
isVal = psY(varargin, 'val', 'n');
vals = ps(varargin, 'vals', []);
yLim = ps(varargin, 'yLim', []);
algs = ps(varargin, 'algs', []);
xs = ps(varargin, 'xs', []);

ticks = ps(varargin, 'ticks', []);
tickAx = ps(varargin, 'tickAx', []);

% body
hold on;

% dimension
[k, n] = size(Me);

% mk & cl
[~, cls] = genMkCl;

% x position
if isempty(xs)
    xs = 1 : n;
end
gap = xs(2) - xs(1);

% width
barWidG2 = barWidG / 2;
barWid = barWid * (barWidG / k);
gapWid = (barWidG - barWid * k) / (k - 1);
devWid2 = barWid * devWid / 2;

% mean
for i = 1 : n
    x = xs(i) - barWidG2;
    for c = 1 : k
        xL = x + (gapWid + barWid) * (c - 1);
        xR = xL + barWid;
        y = Me(c, i);

        miy = 1e-4;
        if bdWid == 0
            fill([xL xR xR xL], [miy miy y y], cls{c}, 'EdgeColor', cls{c});
        else
            fill([xL xR xR xL], [miy miy y y], cls{c}, 'EdgeColor', 'w', 'LineWidth', bdWid);
        end
    end
end
may = max(max(Me));

% legend
if ~isempty(algs)
    h.leg = legend(algs{:}, 'orientation', 'vertical');
end

% standard deviation
if isDev
    for i = 1 : n
        x = xs(i) - barWidG2;

        for c = 1 : k
            xL = x + (gapWid + barWid) * (c - 1);
            xR = xL + barWid;
            xCen = xL + barWid / 2;
            
            y = Me(c, i);
            dev = Dev(c, i);

            xL = xCen - devWid2;
            xR = xCen + devWid2;

            if y == 0
                continue;
            end
        
            % color
            cl = cls{c};

            % vertical line
            plot([xCen xCen], [0 dev] + y, 'Color', cl, 'LineWidth', lnWid);

            % horizontal line
            plot([xL xR], [dev dev] + y, 'Color', cl, 'LineWidth', lnWid); 
        end
    end
end

% boundary
xlim([xs(1) - 1, xs(end) + 1]);
if isempty(yLim)
    mi = 0;
    ma = may * 1.1;
else
    mi = yLim(1); ma = yLim(2);
end
if abs(mi - ma) > eps
    ylim([mi, ma]);
end

% value
if isVal
    for i = 1 : n
        text('Position', [xs(i), max(Me(:, i)) + may * .1], 'String', vals{i}, ...
             'HorizontalAlignment', 'center', 'FontSize', 7);
    end
end

% store
ha.xs = xs;
ha.may = may;

% if ~isempty(ticks)
%     set(gca, 'XTick', []);
%     set(gcf, 'CurrentAxes', tickAx);
%     for i = 1 : n
%         wei = select(true, 'bold', 'normal');
%         text('Position', [i / (n + 1), .5], 'String', ticks{i}, 'Color', 'k', ...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', wei);
%     end
% end
