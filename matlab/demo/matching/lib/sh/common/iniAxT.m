function [axs, axTs, axRs, axTis, axTTis, axRTis] = iniAxT(fig, rows, cols, siz, varargin)
% Create the axes of one figure with title.
%
% Input
%   fig       -  figure number
%   rows      -  number of rows
%   cols      -  number of columns
%   siz       -  figure size, [w, h]
%   varargin
%     pos     -  position of the axes in the figure, {[0 0 1 1]}
%     wGap    -  gap in width, {.2}
%     hGap    -  gap in height, {.2}
%     hT      -  height of top window, {0}
%     wR      -  width of right window, {0}
%     hTi     -  height of main tick window, {0}
%     hTTi    -  height of top tick window, {0}
%     hRTi    -  height of right tick window, {0}
%     ax      -  flag of showing axis, {'y'} | 'n'
%     axT     -  flag of showing axis for text, 'y' | {'n'}
%     name    -  figure name, {''}
%     clf     -  clf flag, {'y'} | 'n'
%     bkCl    -  background color, {'w'}
%
% Output
%   axs       -  handle set of the main window, rows x cols (cell)
%   axTs      -  handle set of the top window, rows x cols (cell)
%   axRs      -  handle set of the right window, rows x cols (cell)
%   axTis     -  handle set of the main tick window, rows x cols (cell)
%   axTTis    -  handle set of the top tick window, rows x cols (cell)
%   axRTis    -  handle set of the right tick window, rows x cols (cell)
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 01-01-2009
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
pos = ps(varargin, 'pos', [0 0 1 1]);
wGap = ps(varargin, 'wGap', .2);
hGap = ps(varargin, 'hGap', .2);
hT = ps(varargin, 'hT', 0);
wR = ps(varargin, 'wR', 0);
hTi = ps(varargin, 'hTi', 0);
hTTi = ps(varargin, 'hTTi', 0);
hRTi = ps(varargin, 'hRTi', 0);
isAx = psY(varargin, 'ax', 'y');
isAxT = psY(varargin, 'axT', 'n');
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
        pos0(3) = siz(1);
        pos0(4) = siz(2);
    end
    
    set(gcf, 'Position', pos0);
end

% sub-figure size
wMar = pos(3) * wGap / (cols + 1); wBody = pos(3) * (1 - wGap) / cols;
hMar = pos(4) * hGap / (rows + 1); hBody = pos(4) * (1 - hGap) / rows;

% locate each axes
[axs, axTs, axRs, axTis, axTTis, axRTis] = cellss(rows, cols);
for row = 1 : rows
    for col = 1 : cols
        % main window
        axs{row, col} = axes('Position', ...
                     [wMar * col + wBody * (col - 1) + pos(1), ...
                      hMar * (rows - row + 1) + hBody * (rows - row + (1 - hT) * hTi) + pos(2), ...
                      wBody * (1 - wR), hBody * (1 - hT) * (1 - hTi)]);
        setAxis(isAx);
        
        % main tick window
        if hTi > 0
            axTis{row, col} = axes('Position', ...
                     [wMar * col + wBody * (col - 1) + pos(1), ...
                      hMar * (rows - row + 1) + hBody * (rows - row) + pos(2), ...
                      wBody * (1 - wR), hBody * (1 - hT) * hTi]);
            setAxis(false);
        end
        
        % top window
        if hT > 0
            axTs{row, col} = axes('Position', ...
                     [wMar * col + wBody * (col - 1) + pos(1), ...
                      hMar * (rows - row + 1) + hBody * (rows - row + (1 - hT) + hT * hTTi) + pos(2), ...
                      wBody * (1 - wR), hBody * hT * (1 - hTTi)]);
            setAxis(isAxT);
        end
        
        % top tick window
        if hT > 0 && hTTi > 0
            axTTis{row, col} = axes('Position', ...
                     [wMar * col + wBody * (col - 1) + pos(1), ...
                      hMar * (rows - row + 1) + hBody * (rows - row + 1 - hT) + pos(2), ...
                      wBody * (1 - wR), hBody * hT * hTTi]);
            setAxis(false);
        end
        
        % right window
        if wR > 0
            axRs{row, col} = axes('Position', ...
                     [wMar * col + wBody * (col - wR) + pos(1), ...
                      hMar * (rows - row + 1) + hBody * (rows - row + (1 - hT) * hRTi) + pos(2), ...
                      wBody * wR, hBody * (1 - hT) * (1 - hRTi)]);
            setAxis(isAxT);  
        end
        
        % right tick window
        if wR > 0 && hRTi > 0
           axRTis{row, col} = axes('Position', ...
                     [wMar * col + wBody * (col - wR) + pos(1), ...
                      hMar * (rows - row + 1) + hBody * (rows - row) + pos(2), ...
                      wBody * wR, hBody * (1 - hT) * hRTi]);
           setAxis(false);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%
function setAxis(flag)
% Set axis.
%
% Input
%   flag  -  axis flag

if flag
    axis on;
else
    axis off;
end 
