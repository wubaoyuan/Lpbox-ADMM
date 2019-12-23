function W=spdiag(D);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

D=D(:);
n=length(D);
W = spdiags(D,0,n,n);


