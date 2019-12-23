function T = tcat(dim, Ts)
% Tensor concatenation.
%
% Input
%   dim     -  dimension
%   Ts      -  tensor set, 1 x m (cell), d_1 x d_2 x d_3
%
% Output
%   T       -  new tensor
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-27-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% original dimension
m = length(Ts);
siz0 = size(Ts{1});
nD = length(siz0);

if nD ~= 3
    error('unsupported');
end

% original specific dimension
ds = zeros(1, m);
for i = 1 : m
    ds(i) = size(Ts{i}, dim);
end
s = n2s(ds);

% new dimension
siz = siz0;
siz(dim) = sum(ds);

% tensor
T = zeros(siz);

for i = 1 : m
    idx = s(i) : s(i + 1) - 1;
    if dim == 1
        T(idx, :, :) = Ts{i};
    elseif dim == 2
        T(:, idx, :) = Ts{i};
    elseif dim == 3
        T(:, :, idx) = Ts{i};
    else
        error('unsupported');
    end
end

