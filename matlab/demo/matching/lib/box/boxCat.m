function [boxs, boxG, Lns, Sca] = boxCat(ori, box0s, varargin)
% Box concatenation.
% 
% Input
%   ori     -  orientation, 'vert' | 'horz'
%   box0s   -  original box set, 1 x m (cell), 2 x 2
%   varargin
%     gap   -  gap size, {0} | ...
%     siz   -  predefined size, {[]} | 2 x 1
%
% Output
%   boxs    -  new box set, 1 x m (cell), 2 x 2
%   boxG    -  global bounding box, 2 x 2
%   Lns     -  breaking line, 1 x (m - 1) (cell), 2 x 2
%   Sca     -  scale (= siz / siz0), 2 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-02-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-08-2012

% function option
gap = ps(varargin, 'gap', 0);
siz = ps(varargin, 'siz', []);

% dimension
m = length(box0s);

% size
Siz0 = boxSiz(box0s);
if isempty(siz)
    siz = max(Siz0, [], 2);
    Sca = ones(2, m);
else
    Sca = repmat(siz, [1, m]) ./ Siz0;
end

% concatenate
mis = cell(1, m);
Lns = cell(1, m - 1);

% vertically
if strcmp(ori, 'vert')
    mis{1} = [0; 0];

    for i = 2 : m
        mis{i} = mis{i - 1} - [0; (1 + gap) * siz(2)];
    end

    % breaking line
    for i = 1 : m - 1
        pt = mis{i} - [0; (1 + .5 * gap) * siz(2)];
        Lns{i} = [pt - [.2 * siz(1); 0], pt + [1.2 * siz(2); 0]];
    end

    boxG = [mis{1} + siz, mis{m}];

% horizontally    
elseif strcmp(ori, 'horz')

    % width & height
    ws = Siz0(1, :);
    hs = Siz0(2, :);
    w = max(ws);
    h = max(hs);
    wGap = w * gap;
    wGap2 = wGap / 2;
    
    % the left-lower corner of each image
    mis{1} = [0; 0];
    for i = 2 : m
        mis{i} = mis{i - 1} + [ws(i - 1) + wGap; 0];
    end
    for i = 1 : m
        mis{i}(2) = (h - hs(i)) / 2;
    end

    % breaking line
    for i = 1 : m - 1
        pt = mis{i} + [ws(i) + wGap2; 0];
        Lns{i} = [pt - [0; .2 * h], pt + [0; 1.2 * h]];
    end

    boxG = [[0; 0], mis{m} + siz];
    
else
    error('unknown ori value: %s', ori);
end

% new box
boxs = cell(1, m);
for i = 1 : m
    boxs{i} = [mis{i}, mis{i} + Siz0(:, i) .* Sca(:, i)];
end
