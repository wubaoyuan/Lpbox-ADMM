function Ax = iniAx(fig, rows, cols, siz, varargin)
% Create axes in a figure.
%
% Input
%   fig     -  figure number
%   rows    -  #rows
%   cols    -  #columns
%   siz     -  figure size, [height, width]
%   varargin
%     pos   -  position of the axes in the figure, {[0 0 1 1]}
%     wGap  -  gap in width, {.2}
%     hGap  -  gap in height, {.2}
%     ws    -  size in width, {[]}
%     hs    -  size in height, {[]}
%     ax    -  flag of showing axis, {'y'} | 'n'
%     name  -  figure name, {''}
%     clf   -  clf flag, {'y'} | 'n'
%     bkCl  -  background color, {'w'}
%
% Output
%   Ax      -  handle set of the main window, rows x cols (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-01-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-04-2012

% function option
pos = ps(varargin, 'pos', [0 0 1 1]);
wGap = ps(varargin, 'wGap', .2);
hGap = ps(varargin, 'hGap', .2);
ws = ps(varargin, 'ws', []);
hs = ps(varargin, 'hs', []);
isAx = psY(varargin, 'ax', 'y');
name = ps(varargin, 'name', '');
isClf = psY(varargin, 'clf', 'y');
bkCl = ps(varargin, 'bkCl', 'w');

% figure
if fig > 0
    figure(fig);

    if isClf
       clf('reset');
    end
    
    set(gcf, 'Color', bkCl);
end

% figure name
if ~isempty(name)
    set(gcf, 'name', name);
end

% figure position & size
if ~isempty(siz)
    if length(siz) == 4
        pos0 = siz;
    else
        pos0 = get(gcf, 'Position');
        pos0(3) = siz(2);
        pos0(4) = siz(1);
    end

    % on Windows
    if ispc
        pos0(1 : 2) = [100 100];
    end
    
    set(gcf, 'Position', pos0);    
end

% axes size
wMar = pos(3) * wGap / (cols + 1);
if isempty(ws)
    wBodys = zeros(1, cols) + pos(3) * (1 - wGap) / cols;
else
    wBodys = ws ./ sum(ws) * (1 - wGap) * pos(3);
end
wws = cumsum([0 wBodys]);

hMar = pos(4) * hGap / (rows + 1);
if isempty(hs)
    hBodys = zeros(1, rows) + pos(4) * (1 - hGap) / rows;
else
    hBodys = hs ./ sum(hs) * (1 - hGap) * pos(4);
end
hhs = cumsum([0 hBodys(end : -1 : 1)]);

% create axes
Ax = cell(rows, cols);
for row = 1 : rows
    for col = 1 : cols
        Ax{row, col} = axes('Position', ...
                     [wMar * col + wws(col) + pos(1), ...
                      hMar * (rows - row + 1) + hhs(rows - row + 1) + pos(2), ...
                      wBodys(col), hBodys(row)]);
        if ~isAx
            set(gca, 'visible', 'off');
        end
    end
end
