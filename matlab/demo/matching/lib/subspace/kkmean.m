function [G, Gs] = kkmean(K, G0, varargin)
% Kernel k-means.
%
% Input
%   K       -  kernel matrix, [] | n x n
%              if K == [], then a KG has been specified in the global environment
%   G0      -  initial indicator matrix, k x n
%   varargin
%     itMa  -  {100}
%
% Output
%   G       -  class indicator matrix, k x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-03-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% kernel in the global environment
global KG;
isKG = isempty(K);

% function option
itMa = ps(varargin, 'itMa', 100);

% dimension
[k, n] = size(G0);
D = zeros(k, n);

Gs = cell(1, itMa + 1);
Gs{1} = G0;
for it = 1 : itMa
    
    % sample-class distance
    for c = 1 : k
        pc = find(Gs{it}(c, :)); 
        nc = length(pc);

        if isKG
            a = sum(sum(KG(pc, pc))) / nc ^ 2;
            for i = 1 : n
                D(c, i) = KG(i, i) - 2 * sum(KG(i, pc)) / nc + a;
            end
        else
            a = sum(sum(K(pc, pc))) / nc ^ 2;
            for i = 1 : n
                D(c, i) = K(i, i) - 2 * sum(K(i, pc)) / nc + a;
            end
        end
    end

    % assign
    [d, l] = min(D);
    
%     Gs{it + 1} = L2G(l, k);
    Gs{it + 1} = zeros(k, n);
    for i = 1 : length(l)
        Gs{it + 1}(l(i), i) = 1;
    end

    % stop condition
%     if cluEmp(Gs{it + 1}) || isequal(Gs{it + 1}, Gs{it})
    if isequal(Gs{it + 1}, Gs{it})
        break;
    end
end
Gs(it + 1 : end) = [];

G = Gs{end};
