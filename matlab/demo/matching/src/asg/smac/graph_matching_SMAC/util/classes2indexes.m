function [indexes,L]=classes2indexes(classes,k);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if nargin<2
    k=max(classes(:));
end
[indexes,L] = mex_classes2indexes(classes,k);
%{
[p,q]=size(classes);
n=p*q;
k=max(classes(:));

ind=find(classes(:));
temp=sparse(1:length(ind),classes(ind),1,n,k);
for i=k:-1:1
    indexes{i,1}=find(temp(:,i));
end

if nargout>=2
    L=sum(temp,1);
    L=full(L(:));
end
%}
