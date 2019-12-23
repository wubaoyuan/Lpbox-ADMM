function saveFigure(figureNumber,isImage,filename);
% TODO : handle function in case remote server without graphics => option
% to save to a given path,filename (works with saveas or print)
% TODO : tiff option to create a preview does not seem to work
% zbuffer, r300, axis image, -deps2, -djpeg
% use getfrme
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.


if nargin < 1
    figureNumber = gcf;
end

if nargin < 2
    isImage=0;
end

if nargin >= 3
    if isImage
        saveFigure_image(figureNumber,filename);
    else
        saveFigureAux(figureNumber,filename);
    end
    return;
end


global data;
if ~isfield(data,'paths')
    data.paths.repertoireResults=pwd;
end
if ~isfield(data,'fileExtension')
    data.fileExtension='.jpg';
end


filepath = fullfile(data.paths.repertoireResults,[get(figureNumber,'Name') data.fileExtension]);
filepath=uiputfile2(filepath,'Save current figure');


[path,name,ext]=fileparts(filepath);
if isempty(ext)
    filepath=fullfile(path,[name,data.fileExtension]);
end

if ~isempty(filepath)
    if isImage
        saveFigure_image(figureNumber,filepath);
    else
        saveFigureAux(figureNumber,filepath);
    end
    [data.paths.repertoireResults,fileName,data.fileExtension] = fileparts(filepath);
    hyperlink(filepath);    
end

function saveFigure_image(fig,filepath);
image=getimage(fig);
% image=getimage(imgca);
cmap=get(gcf,'Colormap');
saveImage(filepath,image,cmap);

