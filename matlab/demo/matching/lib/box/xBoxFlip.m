function X = xBoxFlip(X0, box)
% Flip X.
% 
% Input
%   X0       -  original position, 2 x n | 2 x m x n
%   box      -  new box, 2 x 2
%
% Output
%   X        -  new position, 2 x n | 2 x m x n
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 08-31-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
ds = size(X0);
ndim = length(ds);

% new position
X = X0;

if ndim == 2
    X(2, :) = box(2, 2) - X0(2, :);
elseif ndim == 3
    X(2, :, :) = box(2, 2) - X0(2, :, :);
else
    error('unsupported ndim: %d', ndim);
end
