function ind = sub2ind2(pq,x,y);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if nargin<3
    if isempty(x)
        ind=x;
        return;
    end

    y=x(:,2);
    x=x(:,1);
end
ind = x+pq(1)*(y-1);