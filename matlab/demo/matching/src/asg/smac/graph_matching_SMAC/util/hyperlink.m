function string=hyperlink(path);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.


% disp(['<a href="matlab: winopen(','''',path,'''','); ">',path,'</a>'])

[path0,name,ext]=fileparts(path);
if strcmp(path0,path)
    isDir=1;
else
    isDir=0;
end


if ~isDir
    link1=['<a href="matlab: open2(','''',path,'''','); ">',path,'</a>'];
    link2=['<a href="matlab: open2(','''',path0,'''','); ">','dir','</a>'];
    string=[link1,' ; ',link2];    
else
    string=['<a href="matlab: open2(','''',path,'''','); ">',path,'</a>'];
end


if nargout==0
    disp(string);
end

