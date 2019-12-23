function path = getPaths(name)
% global data;
% path =  getfield(data.paths,name);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.



path='';

if 1%isempty(data) || ~isfield(data,'paths')
    persistent data2;
    if isempty(data2) || nargin<1 %allows to clear
        [data2.paths,data2.user,data2.host] = definePaths();
    end
    if nargin>=1
        path =  data2.paths.(name);
    end
else
    path =  data.paths.(name);
end


