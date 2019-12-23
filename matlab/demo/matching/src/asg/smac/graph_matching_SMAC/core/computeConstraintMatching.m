function [C,b] = computeConstraintMatching(W12, constraint);
% constraint = 'col' | 'row' | {'both'}
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.


if nargin < 2
    constraint = 'both';
end

[n1, n2] = size(W12);
n12 = nnz(W12);
indj = 1 : n12;
indi = zeros(n12, 1);
temp = sum(W12, 1);
starts = cumsum(temp) - temp + 1;
ends = cumsum(temp);
for j = 1 : n2
    indi(starts(j) : ends(j)) = j;
end

% C1 = sparse(indi,indj,1,n2+n1,n12);
C1 = sparse(indi, indj, 1, n2, n12);

indj = 1 : n12;
indi = zeros(n12, 1);
for j = 1 : n2
    indi(starts(j) : ends(j)) = find(W12(:, j));
end

% C2 = sparse(indi+n2,indj,1,n2+n1,n12);
C2 = sparse(indi,indj,1,n1,n12);

% b1 = ones(n2, 1) / n2;
% b2 = ones(n1, 1) / n1;
b1 = ones(n2, 1) / n2 * sqrt(n1 * n2);
b2 = ones(n1, 1) / n1 * sqrt(n1 * n2);

switch constraint
    case 'col'
        C=C1;
        b=b1;
    case 'row'
        C=C2;
        b=b2;
    case 'both'
        C = [C1;C2];
        b = [b1;b2];
end
