function pts1 = ctps_warp_pts(pts, x, c, d)
% Given TPS [t,c,d], warp pts x --> x1.
%           (z is the basis point set).
%           (x1 is in normal coordinate, not expanded !).
%
% Usage: 
%   pts1 = ctps_warp_pts(pts, x, c, d);
%   pts1 = ctps_warp_pts(pts, x, y, lamda);
%   pts1 = ctps_warp_pts(pts, x, y); lamda default = 1;
%
% Input
%   pts     -
%   x       -
%   c       -
%   d       -
%
% Output
%   pts1    -
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-20-2012

% dimension
[dim1, dim2] = size(d); 
% d     -- affine: dim > 1
% lamda -- par:    dim = 1.

% pts1 = ctps_warp_pts(pts, x, c, d) --------------------------
if dim1 > 1
    K = ctps_gen(pts, x);

    % warp pts --> pts1
    [n, dim] = size(pts);
    pts = [ones(n, 1), pts];
    pts1 = pts * d + K * c;
    pts1 = pts1(:, 2 : dim + 1);

% pts1 = ctps_warp_pts(pts, x, y, lamda) ---------------------------
else
    x = x;
    y = c; % y is taking c's position as input now.

    if (nargin == 3)
        lamda = 1;
    else
        lamda = d; % lamda is taking d's position as input now.
    end
    
    K = ctps_gen(pts, x);
    [c, d] = ctps_gen(x, y, lamda);

    % warp pts --> pts1:
    [n, dim] = size(pts);
    pts = [ones(n, 1), pts];
    pts1 = pts * d + K * c;
    pts1 = pts1(:, 2 : dim + 1);
end
