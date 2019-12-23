function target=ind2perm(perm);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

n=length(perm);
if max(perm(:))~=n
    error('TODO');
end
n1=n;
n2=n;
target=full(sparse(1:n1,perm,1,n1,n2))';