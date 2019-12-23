function vs = imgPtVal(R, Pt)
% Compute the response map on the specified position.
%
% Input
%   R       -  response map, h x w
%   Pt      -  points, 2 x k
%
% Output
%   vs      -  scores, k x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-10-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-26-2012

% dimension
[h, w] = size(R);
k = size(Pt, 2);

% float -> integer
Pt2 = round(Pt);
xs = Pt2(1, :);
ys = Pt2(2, :);

% index
idx = sub2ind([h w], ys, xs);

% fetch the value
vs = R(idx);
