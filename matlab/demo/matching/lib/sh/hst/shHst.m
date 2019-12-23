function ha = shHst(mes, devs, varargin)
% Plot histogram in 2-D.
%
% Input
%   mes       -  histogram, 1 x n
%   devs      -  deviation, 1 x n
%   varargin
%     show option
%     barWid  -  width of histogram bar, {.8} <= 1
%     devWid  -  width of deviation line in terms of the 'barWid', {1} <= 1
%     bdWid   -  width of boundary of bar, {1}
%     lnWid   -  width of deviation line, {2}
%     dev     -  dev flag, {'y'} | 'n'
%     val     -  value flag, 'y' | {'n'}
%     form    -  value form, {'%d'}
%     yLim    -  threshold in y axis, {[]}
%     leg     -  legend, {[]}
%
% Output
%   ha        -  handle
%
% Hstory
%   create    -  Feng Zhou (zhfe99@gmail.com), 03-03-2009
%   modify    -  Feng Zhou (zhfe99@gmail.com), 04-17-2013

% show option
psSh(varargin);

% function option
barWid = ps(varargin, 'barWid', .8); 
devWid = ps(varargin, 'devWid', .6);  
bdWid = ps(varargin, 'bdWid', 1);
lnWid = ps(varargin, 'lnWid', 2);
isDev = psY(varargin, 'dev', 'y');
isVal = psY(varargin, 'val', 'n');
form = ps(varargin, 'form', '%d');
yLim = ps(varargin, 'yLim', []);
leg = ps(varargin, 'leg', []);
cls = ps(varargin, 'cls', []);
xs = ps(varargin, 'xs', []);

% body
hold on;

% dimension
n = length(mes);

% mk & cl
if isempty(cls)
    [~, cls] = genMkCl;
end

% x position
if isempty(xs)
    xs = 1 : n;
end
gap = xs(2) - xs(1);

% width
barWid2 = barWid / 2 * gap;
devWid2 = devWid / 2 * gap;

% mean
for i = 1 : n
    x = xs(i);
    y = mes(i);    
    xL = x - barWid2;
    xR = x + barWid2;

    if y == 0
        continue;
    end

    % color
    cl = cls{1 + mod(i - 1, length(cls))};

    if bdWid == 0
        fill([xL xR xR xL], [0 0 y y], cl, 'EdgeColor', cl);
    else
        fill([xL xR xR xL], [0 0 y y], cl, 'LineWidth', bdWid);
    end
end
ys = mes;

% legend
if ~isempty(leg)
    h.leg = legend(leg{:}, 'orientation', 'vertical');
end

% standard deviation
if isDev && ~isempty(devs)
    ys = mes + devs;

    for i = 1 : n
        x = xs(i);
        y = mes(i);
        dev = devs(i);
        xL = x - devWid2;
        xR = x + devWid2;

        if y == 0
            continue;
        end
        
        % color
        cl = cls{1 + mod(i - 1, length(cls))};

        % vertical line
        plot([x x], [0 dev] + y, 'Color', cl, 'LineWidth', lnWid);

        % horizontal line
        plot([xL xR], [dev dev] + y, 'Color', cl, 'LineWidth', lnWid); 
    end
end

% maximum value
may = max(ys);

% value
if isVal
    [val, vals] = vec2str(mes, form);
    for c = 1 : k
        text('Position', [xs(c), mes(c) + may * .1], 'String', vals{c}, ...
            'HorizontalAlignment', 'center');
    end
end

% boundary
xlim([xs(1) - gap, xs(end) + gap]);
if isempty(yLim)
    mi = 0;
    ma = may * 1.1;
else
    mi = yLim(1); ma = yLim(2);
end
if abs(mi - ma) > eps
    ylim([mi, ma]);
end

% store
ha.xs = xs;
ha.ys = ys;
