function h=getimageHandle(h0);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if nargin<1
    h0=gca;
end

hs=get(h0,'Children');
for i=1:length(hs)
    hi=hs(i);
    if strcmp(get(hi,'Type'),'image')
        h=hi;
        return;
    end
end

h=[];
