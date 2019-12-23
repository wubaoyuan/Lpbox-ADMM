function s=emptyStruct(fields,sizeStruct);
%{
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

res(nFiles,1)=emptyStruct([]);
res=emptyStruct(fields,[10,10]);
%}

n=length(fields);
s=cell(2*n,1);
s(1:2:end)=fields(:);

s2={{}};
s2=s2(ones(n,1));
s(2:2:end)=s2;

s=struct(s{:});
if nargin>=2
    if isempty(s)
        s(prod(sizeStruct)).(fields{1})=[];
    else
        s(prod(sizeStruct))=s;%doesn't work since s is empty
    end
    s=reshape(s,sizeStruct);
end

