function X = xBoxNew(X0, box0, box, sca)
% Obtain new X.
% 
% Input
%   X0       -  original position, 2 x n | 2 x m x n
%   box0     -  original box, 2 x 2
%   box      -  new box, 2 x 2
%   sca      -  scale (= siz / siz0), 2 x 1
%
% Output
%   X        -  new position, 2 x n | 2 x m x n
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 08-31-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
ds = size(X0);
ds(1) = 1;

% new position
X = (X0 - repmat(box0(:, 1), ds)) .* repmat(sca, ds) + repmat(box(:, 1), ds);
