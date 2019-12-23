function [ns, idx] = divN(n, m, varargin)
% Decompose a set of n points into m subsets according the specified way.
% The size of subset must sum to n. i.e. sum(ns) == n.
%
% Example
%   input   -  n = 10; m = 3; alg = 'eq';
%   call    -  [ns, idx] = divN(n, m, 'alg', alg)
%   output  -  ns = [3; 3; 4];
%              idx = [2; 1; 2; 3; 3; 3; 1; 2; 3; 1];
%
% Input
%   n       -  number
%   m       -  parts
%   varargin
%     alg   -  algorithm type, 'rand' | 'w' | {'eq'}
%
% Output
%   ns      -  size of subset, m x 1
%   idx     -  index of subset that each number (1 ~ n) has been assigned to, n x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 05-30-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
alg = ps(varargin, 'alg', 'eq');

if strcmp(alg, 'eq')
    ns = zeros(m, 1);
    w = floor(n / m);
    for i = 1 : m - 1
        ns(i) = w;
    end
    ns(m) = n - w * (m - 1);

elseif strcmp(alg, 'w')
    w = m;
    m = ceil(n / w);
    ns = zeros(m, 1);
    for i = 1 : m - 1
        ns(i) = w;
    end
    ns(m) = n - w * (m - 1);
    
    if ns(m) < w / 2
        ns(m - 1) = ns(m - 1) + ns(m);
        ns(m) = [];
    end
    
elseif strcmp(alg, 'rand')
    while true
        aa = rand(m, 1); aa = aa / sum(aa);

        ns = floor(n * aa);
        n2 = n - sum(ns);

        i = 0;
        while n2 > 0
            i = mod(i + 1, m) + 1;
            ns(i) = ns(i) + 1;

            n2 = n2 - 1;
            i = i + 1;
        end

        if isempty(find(ns == 0, 1))
            break;
        end
    end
else
    error(['unknown algorithm: ' alg]);
end

% index
idx0 = randperm(n);
idx = zeros(n, 1);
cns = [0; cumsum(ns)];
for c = 1 : m
    idx(idx0(cns(c) + 1 : cns(c + 1))) = c;
end

