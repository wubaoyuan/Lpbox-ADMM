function path2=changeExt(path1,ext2);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if ~strcmp(ext2(1),'.')
    error('ext2 must start with .');
end
[path,name,ext]=fileparts(path1);
path2=fullfile(path,[name,ext2]);

