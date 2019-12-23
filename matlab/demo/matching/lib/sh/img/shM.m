function ha = shM(M, varargin)
% Show matrix in 2-D space.
%
% Input
%   M        -  matrix, n1 x n2
%   varargin
%     show option
%     dis    -  display function, {'imagesc'} | 'contour'
%     clMap  -  color map, {'gray'}
%     bar    -  bar flag, 'y' | {'n'}
%     eq     -  axis equal flag, 'y' | {'n'}
%     P      -  warping path, {[]}
%     lnWid  -  line width, {1}
%     lnCl   -  line color (for boundary), {'r'}
%
% Output
%   ha       -  figure handle
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify   -  Feng Zhou (zhfe99@gmail.com), 04-24-2013

% show option
psSh(varargin);

% function option
dis = ps(varargin, 'dis', 'imagesc');
clMap = ps(varargin, 'clMap', 'gray');
isBar = psY(varargin, 'bar', 'n');
isEq = psY(varargin, 'eq', 'y');
isSq = psY(varargin, 'sq', 'n');
P = ps(varargin, 'P', []);
Ps = ps(varargin, 'Ps', []);
lnMk = ps(varargin, 'lnMk', '-');
lnWid = ps(varargin, 'lnWid', 1);
lnCl = ps(varargin, 'lnCl', 'r');
cmap = ps(varargin, 'cmap', []);
clim = ps(varargin, 'clim', []);

% matrix
[n1, n2] = size(M);
if strcmp(dis, 'imagesc')
    if strcmp(clMap, 'grayc')
        M = repmat(M, [1 1 3]);
    end
    if isempty(clim)
        ha.M = imagesc(M);
    else
        ha.M = imagesc(M, clim);
    end
elseif strcmp(dis, 'contour')
    ha.M = contour(M);
else
    error(['unknown display type: ' dis]);
end

% color map
if ~isempty(cmap)
    ha.cmap = colormap(cmap);
elseif strcmp(clMap, 'gray')
    ha.cmap = colormap(gray);
elseif strcmp(clMap, 'grayc')
    ha.cmap = colormap('jet');
elseif strcmp(clMap, 'hsv')
    ha.cmap = colormap(hsv);
elseif strcmp(clMap, 'jet')
    ha.cmap = colormap(jet);
elseif strcmp(clMap, 'hot')
    ha.cmap = colormap(hot);    
else
    error('unknown color map: %s', clMap);
end

% axis
if isEq
    axis equal;
end
if isSq
    axis square;
end
axis([1 - .5, n2 + .5, 1 - .5, n1 + .5]);
axis ij;
hold on;
set(gca, 'ticklength', [0 0]);
%set(gca, 'XTick', 1 : n2, 'YTick', 1 : n1);

% color bar
if isBar
    colorbar;
end

% warping path
if ~isempty(P)
    plot(P(:, 2), P(:, 1), lnMk, 'Color', lnCl, 'LineWidth', lnWid);
end
if ~isempty(Ps)
    for i = 1 : length(Ps)
        P = Ps{i};
        plot(P(:, 2), P(:, 1), '--', 'Color', 'g', 'LineWidth', lnWid);
    end    
end
