function [C, b] = compute_full_rank_constraint(C, b);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

[k, n] = size(C);
% assert(k<n);
r = sprank(C);

%TODO:caution:sprank(A)>=rank(full(A))... adapt?

if r < k
    [b, C] = qr(C, b);
    temp = sum(C .^ 2, 2);
    ind = find(temp > eps);
    C = C(ind, :);
    b = b(ind);
end
