function vs = imgPeakBox(Rs, Pt)
% Compute the response map on the specified position.
%
% Input
%   Rs      -  response map, h x w x k
%   Pt      -  points, 2 x k
%
% Output
%   vs      -  scores, k x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-10-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-16-2012

% dimension
[h, w, k] = size(Rs);

% float -> integer
Pt2 = round(Pt);

% valid position
iXs = Pt2(1, :);
iYs = Pt2(2, :);
vis = iXs >= 1 & iXs <= w & iYs >= 1 & iYs <= h;

% index
cs = 1 : k;
idx = sub2ind([h w k], iYs(vis), iXs(vis), cs(vis));

% fetch the value
vs = zeros(k, 1);
vs(vis) = Rs(idx);
