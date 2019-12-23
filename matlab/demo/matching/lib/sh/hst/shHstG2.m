function h = shHstG(HistG, varargin)
% Plot histogram group in 2-D.
%
% Input
%   HistG      -  histogram group matrix, k x kG x n
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
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 03-03-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% show option
psShow(varargin);

% function option
barWidG = ps(varargin, 'barWidG', .8); 
barWid = ps(varargin, 'barWid', .8);    
devWid = ps(varargin, 'devWid', .8);    
bdWid = ps(varargin, 'bdWid', 1);
lnWid = ps(varargin, 'lnWid', 2);
isDev = psY(varargin, 'dev', 'y');
isVal = psY(varargin, 'val', 'n');
form = ps(varargin, 'form', '%d');
yLim = ps(varargin, 'yLim', []);
leg = ps(varargin, 'leg', []);
ticks = ps(varargin, 'ticks', []);
tickAx = ps(varargin, 'tickAx', []);

% body
hold on;

% class
[markers, colors] = genMarkers;
[k, kG, n] = size(HistG);

% width
barWidG2 = barWidG / 2;
barWid = barWid * (barWidG / k);
gapWid = (barWidG - barWid * k) / (k - 1);
devWid = devWid * barWid;
xs = 1 : kG;

% mean
Me = mean(HistG, 3);
for cG = 1 : kG
    x = xs(cG) - barWidG2;
    for c = 1 : k
        xL = x + (gapWid + barWid) * (c - 1);
        xR = xL + barWid;
        y = Me(c, cG);

        if bdWid == 0
            fill([xL xR xR xL], [0 0 y y], colors{c}, 'EdgeColor', colors{c});
        else
            fill([xL xR xR xL], [0 0 y y], colors{c}, 'EdgeColor', 'w', 'LineWidth', bdWid);
        end
    end
end
may = max(max(Me));

% legend
if ~isempty(leg)
    h.leg = legend(leg{:}, 'orientation', 'vertical');
end

% standard deviation
if isDev
    Dev = std(HistG, 0, 3);
    for cG = 1 : kG
        x = xs(cG) - barWidG2;
        for c = 1 : k
            xM = x + (gapWid + barWid) * (c - 1) + barWid / 2;
            xL = xM - devWid / 2;
            xR = xM + devWid / 2;
            y = Me(c, cG);
            dev = Dev(c, cG);
            
            cl = colors{c};

            % vertical line
            line([xM xM], [0 dev] + y, 'Color', cl, 'LineWidth', lnWid);

            % horizontal line
            line([xL xR], [dev dev] + y, 'Color', cl, 'LineWidth', lnWid); 
        end
    end
    may = may + max(max(Dev));
end

disp(Me);
if isDev
    disp(Dev);
end

% value
if isVal
%     [val, vals] = vec2str(mes, form);
%     for cG = 1 : kG
%         for c = 1 : k
%             text('Position', [xs(c), mes(c) + may * .1], 'String', vals{c}, ...
%                 'HorizontalAlignment', 'center');
%         end
%     end
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

if ~isempty(ticks)
    set(gca, 'XTick', []);
    set(gcf, 'CurrentAxes', tickAx);
    for i = 1 : n
        wei = select(true, 'bold', 'normal');
        text('Position', [i / (n + 1), .5], 'String', ticks{i}, 'Color', 'k', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', wei);
    end
end
