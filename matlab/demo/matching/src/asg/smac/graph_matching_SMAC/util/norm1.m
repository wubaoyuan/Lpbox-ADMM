function z = norm1(X);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

X=double(X);
if prod(size(X)) < 2^31-2
    z = norm(X(:));
else
    z = sqrt(sum(vec(sum(abs(X).^2,2))));
end

z=full(z);
