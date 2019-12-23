function h = show3dIni(X, H, varargin)
% Show data as points with different labels for different classes in 3-D figure.
%
% Input
%   X        -  sample matrix, dim x n
%   H        -  indicator matrix, k x n. If empty, H is considered as ones(1, n)
%   varargin
%     mkSiz  -  size of markers, {5}
%     leg    -  legends, {[]}
%     anotI  -  index of point, {[]}
%     anotS  -  string of annotation, {[]}
%     marks  -  markers, {[]}
%     face   -  flag of showing face color, 'y' | {'n'}
%
% Output
%   h        -  figure content handle
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 12-30-2008
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
psShow(varargin);
mkSiz = ps(varargin, 'mkSiz', 5);
leg = ps(varargin, 'leg', []);
anotI = ps(varargin, 'anotI', []);
anotS = ps(varargin, 'anotS', []);
marks = ps(varargin, 'marks', []);
isFace = psY(varargin, 'face', 'n');

% dimension
[dim, n] = size(X);
if dim ~= 3
    error('incorrect dim');
end

% label
if isempty(H), H = ones(1, n); end
k = size(H, 1); L = G2L(H);

% markers
[markers, colors] = genMarkers;
if k > length(markers)
    error('too much classes');
end
if ~isempty(marks)
    markers = marks;
end

% point
hold on;
poi = cell(1, k);
for c = 1 : k
    Y = X(:, L == c);
    
    poi{c} = plot3(Y(1, :), Y(2, :), Y(3, :), markers{c}, ...
        'MarkerSize', mkSiz, 'Color', colors{c}, 'LineWidth', 1);

    if isFace
        set(poi{c}, 'MarkerFaceColor', colors{c})
    end
    
    if ~isempty(leg)
        set(poi{c}, 'DisplayName', leg{c});
    end
end
h.poi = poi;

% legend
h.leg = [];
if ~isempty(leg)
    h.leg = legend(leg{:});
end

% annotation
nA = length(anotI);
if nA > 0 && isempty(anotS)
    [anotTmp, anotS] = vec2str(1 : nA, '%d');
end

for iA = 1 : nA
    ind = anotI(iA);
%     plot(X(1, ind), X(2, ind), 'o', 'Color', 'k', 'MarkerSize', mkSiz + 2);
    text(X(1, ind), X(2, ind), anotS{iA});
end

% boundary
setBound(X, 'mar', [.1 .1 .1]);
