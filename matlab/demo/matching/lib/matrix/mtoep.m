function X = mtoep(alg, x, m)
% Obtain a Toeplitz-structured matrix.
%
% Input
%   alg     -  algorithm type, 'u' | 'a'
%              'u': upper Toeplitz matrix
%              'a': asymmetric Toeplitz matrix
%   x       -  signal, (2n - 1) x 1 | n x 1
%   m       -  parameter, used only when alg == 2
% 
% Output
%   X       -  Toeplitz matrix, n x n | n x (n + m - 1)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-02-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
if strcmp(alg, 'u')
    n = size(x, 1);
    
    X = zeros(m, n + m - 1);

    idx = 1 : m + 1 : (m + 1) * (m - 1) + 1;
    for i = 1 : n
        X(idx) = x(i);
        idx = idx + m;
    end
    
elseif strcmp(alg, 'a')
    n = length(x);

    n = round((n + 1) / 2);

    X = zeros(m, m);

    % bottom-left
    idx = n : m + 1 : (m + 1) * (m - n) + n;
    for i = 1 : n
        X(idx) = x(i);
        idx = [idx - 1, idx(end) + m];
    end

    % upper-right
    idx = 1 : m + 1 : (m + 1) * (m - 1) + 1;
    for i = 2 : n
        idx = idx(1 : end - 1) + m;
        X(idx) = x(i + n - 1);
    end

else
    error('unknown alg: %s', alg);
end
