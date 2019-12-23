function filepath=uiputfile2(filter,query_string);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if nargin<2
    query_string='Choose file name to save to';
end

[fileName,pathName] = uiputfile(filter,query_string);
filepath=[];
if ~(isequal(fileName,0)|isequal(pathName,0) )
    [ignore,fileName,ext] = fileparts(fileName);
    filepath=[pathName fileName ext];    
end

