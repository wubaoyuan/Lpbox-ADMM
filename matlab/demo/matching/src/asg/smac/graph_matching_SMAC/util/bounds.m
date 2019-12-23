function b = bounds(x);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if isempty(x)
    b=[];
    return;
end
if issparse(x)
    b(1)=min(min(x));
    b(2)=max(max(x));
else
    b(1) = min(x(:));
    b(2) = max(x(:));
end
b = full(b);