function D = conDistSpace(Xs, varargin)
% Construct the distance between the subspace defined by the sample matrix.
%
% Input
%   Xs      -  set of sample matrix, 1 x m (cell), dim x n
%   varargin 
%     alg   -  distance type, {'ang'} | 'eig'
%     val   -  dependent on the algorithm, {10} for 'eig'
%
% Output
%   D      -   distance matrix, m x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-30-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
alg = ps(varargin, 'alg', 'ang');
val = ps(varargin, 'val', 10);

m = length(Xs);
D = zeros(m, m);

% svd of each matrix
[Us, Ss] = cellss(1, m);
for i = 1 : m
    [U, S] = svd(Xs{i});
    [Ss{i}, ind] = sort(sum(S, 2), 'descend');
    Us{i} = U(:, ind);
end

for i = 1 : m
    for j = i + 1 : m
        Xi = Xs{i};
        Xj = Xs{j};
        
        if strcmp(alg, 'ang')
            d = subspace(Xi, Xj);

        elseif strcmp(alg, 'eig')
            Si = Ss{i}; Ui = Us{i};
            Sj = Ss{j}; Uj = Us{j};
            
            d = 0;
            for k = 1 : val
                d = d + Si(k) / Sj(k) * Ui(:, k)' * Uj(:, k);
            end

        else
            error('unknown distance type');
        end
        D(i, j) = d;
    end
end
