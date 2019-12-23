function saveImage(filepath,image,cmap);
% saveImage(dataResult,'image',image,gray);
% dataResult.dirOutput
% dataResult.prefix
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

    
% [ignore,name,ext]=fileparts(name);
% if isempty(ext)
%     ext='.png';
% end
[p,q,r]=size(image);
if r==1
    if nargin<3
        image=gray2rgb(image,jet(128));
    else
        image=gray2rgb(image,cmap);
    end
end
% filename=fullfile(dataResult.dirOutput,[dataResult.prefix,'_',name,ext]);
if isempty(image)
    image=0;
end
imwrite(image,filepath);
