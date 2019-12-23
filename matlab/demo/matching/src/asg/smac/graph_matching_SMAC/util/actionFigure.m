function actionFigure(hObject, eventdata, handles);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

global data;
selectionType = get (gcf , 'SelectionType');
fig = gcf;

switch selectionType
    case 'alt'
        cmenu = uicontextmenu;
        fun=@(hObject, eventdata,handles)(add_scroll_panel(fig,hObject, eventdata));
        uimenu(cmenu, 'Label', 'SCROLL','Callback',fun);
        
        fun=@(hObject, eventdata,handles)(saveFigure(fig,0));
        uimenu(cmenu, 'Label', 'saveFigure','Callback',fun);
        fun=@(hObject, eventdata,handles)(saveFigure(fig,1));
        uimenu(cmenu, 'Label', 'saveImage','Callback',fun);
        %         uimenu(cmenu, 'Label', 'open in new window','Callback',@openInNewWindow);
        %         uimenu(cmenu, 'Label', 'analyze data point','Callback',@analysePoint);
        uimenu(cmenu, 'Label', 'axis image','Callback',@axisImage);
        if norm1(get(gca,'Position')-[0,0,1,1])
            uimenu(cmenu, 'Label', 'remove margin','Callback',@removeMargin);
        else
            uimenu(cmenu, 'Label', 'replace margin','Callback',@replaceMargin);
        end

        if strcmp(get(gcf,'MenuBar'),'none')
            fun=@(hObject, eventdata,handles)(set(fig,'MenuBar','figure'));
            uimenu(cmenu, 'Label', 'MenuBar on','Callback',fun);
        else
            fun=@(hObject, eventdata,handles)(set(fig,'MenuBar','none'));
            uimenu(cmenu, 'Label', 'MenuBar off','Callback',fun);
        end
        uimenu(cmenu, 'Label', 'Colorbar','Callback',@colorbar2);
        if norm1(get(gcf,'colormap')-gray)
            uimenu(cmenu, 'Label', 'Colormap gray','Callback',@colormapGray);
        else
            uimenu(cmenu, 'Label', 'Colormap default','Callback',@colormapDefault);
        end
        uimenu(cmenu, 'Label', 'pixel property','Callback',@pixelProperty);
        fun=@(hObject, eventdata,handles)(impixelinfo);
        uimenu(cmenu, 'Label', 'impixelinfo','Callback',fun);
        
        uimenu(cmenu, 'Label', 'getMeasures','Callback','getMeasures()');
        uimenu(cmenu, 'Label', 'afficheGrid','Callback','afficheGrid();');
        uimenu(cmenu, 'Label', 'flip-ud,invert color','Callback',@flipUpDown);
        uimenu(cmenu, 'Label', 'save to data.image','Callback',@saveImageToData);

        if data.gui.dataFigures(fig).isClosable
            fun=@(hObject, eventdata,handles)( eval('data.gui.dataFigures(fig).isClosable=0;') );
            uimenu(cmenu, 'Label', 'set closable=0','Callback',fun);
        else
            fun=@(hObject, eventdata,handles)( eval('data.gui.dataFigures(fig).isClosable=1;') );
            uimenu(cmenu, 'Label', 'set closable=1','Callback',fun);
        end

        UserData=get(gcf,'UserData');
        if isfield(UserData,'h_getPlot')
            uimenu(cmenu, 'Label', 'getPlot','Callback',@getPlot);
        end

        set(gco,'UIContextMenu',cmenu,'Visible','on');%,'Selected','on');
    case 'open'
        openInNewWindow;
    case 'extend'
        % todo : choose what to do
    case 'normal'
        try
            zoomInWindow;
        catch
        end
    otherwise
end


function zoomInWindow(hObject, eventdata,handles)
return;%disabled

global data;
fig0=gcf;
plot0=gca;
if~isfield(data.gui,'zoomWindow')
    data.gui.zoomWindow=figure2;
else
    fig=data.gui.zoomWindow;
    if ~(ishandle(fig) && strcmp(get(fig,'Type'),'figure'))
        data.gui.zoomWindow=figure2;
    end
end
L=data.gui.zoomL;
fig=data.gui.zoomWindow;
if(fig0==fig)
    return;
end
figure(fig);
clf(fig);

set(fig,'name','ZOOM');
set(fig,'MenuBar', 'none');
% set(gcf,'UIContextMenu',cmenu,'Visible','on');
% h=uicontrol(fig,'Style','edit','Position', [0 0 margin 1-margin],'String',num2str(data.gui.zoomL),'Callback',@edit_callback);
h=uicontrol(fig,'Style','edit','String',num2str(data.gui.zoomL),'Callback',@edit_callback);
% UserData.sliderZoom =
% uicontrol(fig,'Style','slider','Units','Normalized',
% set(fig,'UserData',UserData);

CurrentPoint=get(plot0,'CurrentPoint');
x=round(CurrentPoint(1,2));
y=round(CurrentPoint(1,1));
imageW=getimage(plot0);
if isempty(imageW)
    figure(fig0);
    return;
end

[p,q,r]=size(imageW);
newAxes = copyobj(plot0,fig);
set(newAxes,'Units','normalized');
set(newAxes,'Position',[0 0 1 1]);

% subplot(newAxes);
w0=axis(plot0);
w(1)=max(1,x-L);
w(2)=min(p,x+L);
w(3)=max(1,y-L);
w(4)=min(q,y+L);
axis(newAxes,[w(3:4),w(1:2)]);
set(fig,'Colormap',get(fig0,'Colormap'));
figure(fig0);

function zoomInWindow_old(hObject, eventdata,handles)
global data;
fig0=gcf;
plot0=gca;
if~isfield(data.gui,'zoomWindow')
    data.gui.zoomWindow=figure2;
else
    fig=data.gui.zoomWindow;
    if ~(ishandle(fig) && strcmp(get(fig,'Type'),'figure'))
        data.gui.zoomWindow=figure2;
    end
end
L=data.gui.zoomL;
fig=data.gui.zoomWindow;
set(fig,'name','ZOOM');
set(fig,'MenuBar', 'none');
% set(gcf,'UIContextMenu',cmenu,'Visible','on');
% h=uicontrol(fig,'Style','edit','Position', [0 0 margin 1-margin],'String',num2str(data.gui.zoomL),'Callback',@edit_callback);
h=uicontrol(fig,'Style','edit','String',num2str(data.gui.zoomL),'Callback',@edit_callback);
% UserData.sliderZoom =
% uicontrol(fig,'Style','slider','Units','Normalized',
% set(fig,'UserData',UserData);

CurrentPoint=get(plot0,'CurrentPoint');
x=round(CurrentPoint(1,2));
y=round(CurrentPoint(1,1));
imageW=getimage(plot0);
if isempty(imageW)
    figure(fig0);
    return;
end

[p,q,r]=size(imageW);
w(1)=max(1,x-L);
w(2)=min(p,x+L);
w(3)=max(1,y-L);
w(4)=min(q,y+L);
imageW=imageW(w(1):w(2),w(3):w(4),:);
figure(fig);
imagesc(imageW);
axis image
hold on;
plot(min(1+L,y),min(1+L,x),'+k');
plot(min(1+L,y),min(1+L,x),'ow');
hold off;
% set(gca,'Position',[0.1300, 0.1100, 0.7750, 0.8150]);
axis on
figure(fig0);

function edit_callback(hObject, eventdata,handles)
global data;
L=str2num(get(hObject,'String'));
data.gui.zoomL = L;

function openInNewWindow(hObject, eventdata,handles)
copy2figure;

function saveImageToData(hObject, eventdata,handles)
global data;
XLim = round(get(gca,'XLim'));
YLim = round(get(gca,'YLim'));
temp = get(gca,'Children');
temp = get(temp(end));
imageX = double(temp.CData);
[p,q,r] = size(imageX);
i1 = max(1,YLim(1));
i2 = min(p,YLim(2));
j1 = max(1,XLim(1));
j2 = min(q,XLim(2));
imageX = imageX(i1:i2,j1:j2,:);
[data.ImageOriginale , data.imageCouleur , data.image] = computeImagesGui(imageX , data.p,data.q);

%{
function analysePoint(hObject, eventdata,handles)
currentFigure = gcf;

z=get(gca,'CurrentPoint');
UserData = get(gca,'UserData');
X = UserData.X;
Example = UserData.Example;
index = dataPoint2image(Example,X,z(1,:));
figure(currentFigure);
hold on;
selected = num2cell(X(index,:));
if length(selected)<=2
h = plot(selected{:},'or');
else
h = plot3(selected{:},'or');
end
hold off;
%pause(0.2);
%delete(h);
%}

function colorbar2(hObject, eventdata,handles)
%set(gca,'Position',[0 0 0.75 1]);
colorbar;
set(gca,'Position',[0 0 0.75 1]);

function colormapGray(hObject, eventdata,handles)
colormap(gray);

function colormapDefault(hObject, eventdata,handles)
colormap('default');

function axisImage(hObject, eventdata,handles)
axis image

function removeMargin(hObject, eventdata,handles)
set(gca,'Position',[0,0,1,1]);
axis off


function pixelProperty(hObject, eventdata,handles)
pixel=get(gca,'CurrentPoint');
i = round(pixel(1,2));
j = round(pixel(1,1));
image = getimage(gca);
isImage=~isempty(image);
if isImage
    [p,q,r] = size(image);
    if r > 1
        rgb = image(i,j,:);
        hsv(1,1,1:3) = rgb2hsv(rgb);
    else
        intensity = image(i,j);
    end
end
disp(sprintf('\n'));
disp(['[i,j] = [',num2str(i),',',num2str(j),']']);

if isImage
    disp(['z = ',num2str(sub2ind2([p,q],i,j))]);
    if r > 1
        disp(['rgb = [',num2str(rgb(1)),',',num2str(rgb(2)),',',num2str(rgb(3)),']']);
        disp(['hsv = [',num2str(hsv(1)),',',num2str(hsv(2)),',',num2str(hsv(3)),']']);
    else
        disp(['intensity = ',num2str(intensity)]);
    end
    disp(['[p,q,r] = ',num2str([p,q,r])]);
    disp(['bounds(image) = ',num2str(bounds(image))]);
end

h1=dispCircle([i,j],1,'r');
h2=dispCircle([i,j],2,'g');
pause(1);
delete(h1);
delete(h2);

function flipUpDown(hObject, eventdata,handles)
image = getimage(gca);
image = 1-image(end:-1:1,:,:);
imagesc(image);

function getPlot(hObject, eventdata,handles);
UserData=get(gcf,'UserData');
feval(UserData.h_getPlot);

function add_scroll_panel(fig,hObject, eventdata,handles)
h=getimageHandle(gca);
% get(gca,'Children');
hScroll = imscrollpanel(fig,h);
api = iptgetapi(hScroll);

mag = api.findFitMag();
api.setMagnification(mag);

hMagBox = immagbox(fig,h); 
pos = get(hMagBox,'Position');
set(hMagBox,'Position',[0 0 pos(3) pos(4)])
imoverview(h);


