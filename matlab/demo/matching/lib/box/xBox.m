function box = xBox(X, par)
% Obtain bounding box of sample matrix X.
%
% Input
%   X       -  sample matrix, dim x n
%   par     -  function parameter
%     box0  -  pre-specified bounding box, {[]} | dim x 2
%     siz0  -  pre-specified size, {[]} | dim x 1
%     w2h0  -  pre-specified width / height ratio, {[]} | 1 x 1
%     h2w0  -  pre-specified height / width ratio, {[]} | 1 x 1
%     mar   -  pre-specified margin, {.1}
%
% Output
%   box     -  bounding box, dim x 2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 11-20-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function parameter
box0 = ps(par, 'box0', []);
siz0 = ps(par, 'siz0', []);
w2h0 = ps(par, 'w2h0', []);
h2w0 = ps(par, 'h2w0', []);
mar = ps(par, 'mar', .1);

% dimension
dim = size(X, 1);

% bounding box position
if ~isempty(box0)
    box = box0;

% bounding box size
elseif ~isempty(siz0)
    center = mean(X, 2);
    if length(siz0) == 1
        siz = siz0(ones(dim, 1));
    else
        siz = siz0(:);
    end
    box = [center - siz / 2, center + siz / 2];

% w2h ratio
elseif ~isempty(w2h0)
    center = mean(X, 2);
    mi = min(X, [], 2);
    ma = max(X, [], 2);
    siz = ma - mi;
    
    if length(mar) == 1
        mar = ones(dim, 1) * mar;
    else
        mar = mar(:);
    end
    
    siz = siz .* (1 + mar);
    siz2 = siz / 2;
    siz2(2) = w2h * siz2(1);

    box = [center - siz2, center + siz2];

% h2w ratio
elseif ~isempty(h2w0)
    center = mean(X, 2);
    mi = min(X, [], 2);
    ma = max(X, [], 2);
    siz = ma - mi;
    
    if length(mar) == 1
        mar = ones(dim, 1) * mar;
    else
        mar = mar(:);
    end
    
    siz = siz .* (1 + mar);
    siz2 = siz / 2;
    siz2(1) = h2w * siz2(2);

    box = [center - siz2, center + siz2];

% margin
else
    mi = min(X, [], 2);
    ma = max(X, [], 2);
    siz = ma - mi;
    siz(siz == 0) = .1;
    
    if length(mar) == 1
        mar = ones(dim, 1) * mar;
    else
        mar = mar(:);
    end
    
    box = [mi - siz .* mar, ma + siz .* mar];
end
