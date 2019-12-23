function setPropertyFigure2(fig);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

isSlider = 0;

UserData = get(fig,'UserData');

set(fig,'WindowButtonDownFcn',@actionFigure);
% set(fig,'KeyPressFcn',@actionFigure);

if isSlider
    margin = 0.02;
    % UserData.panel = uipanel('Units','Normalized','Position',[0,margin,1-margin,1]);
    % UserData.panel = uipanel('Units','Normalized','Position',[0,-0.3,1.3,11.3]);
    UserData.sliderH = uicontrol(fig,'Style','slider','Units','Normalized', 'Position', [0 0 1-margin margin],'Callback',@slider_callbacks);
    UserData.sliderV = uicontrol(fig,'Style','slider','Units','Normalized', 'Position', [1-margin margin margin 1-margin],'Value',1,'Callback',@slider_callbacks);
    UserData.sliderZoom = uicontrol(fig,'Style','slider','Units','Normalized', 'Position', [0 0 margin 1-margin],'Value',1,'Callback',@slider_callbacks);

    set(fig,'UserData',UserData);
end

function slider_callbacks(hObject, eventdata,handles);
% fig = gcf;
fig = get(hObject,'Parent');
UserData = get(fig,'UserData');
if ~isfield(UserData,'plots')
    UserData.plots = struct('handle',{},'pos0',{});
end

valueZoom = get(UserData.sliderZoom,'Value');
valueZoom = 10^(1-valueZoom);

% set(UserData.sliderH,'SliderStep',[1/valueZoom-1+eps,1/valueZoom-1+eps]);
% set(UserData.sliderV,'SliderStep',[1/valueZoom-1+eps,1/valueZoom-1+eps]);
sliderStep = [1/valueZoom,1/valueZoom];
set(UserData.sliderH,'SliderStep',sliderStep);
set(UserData.sliderV,'SliderStep',sliderStep);
% set(UserData.sliderV,'Max',1);
% set(UserData.sliderH,'Max',1);
valueH = get(UserData.sliderH,'Value');
valueV = get(UserData.sliderV,'Value');
valueV = 1-valueV;

valueH = valueH*valueZoom;
valueV = valueV*valueZoom;
% [valueH,valueV,valueZoom]

children = get(gcf,'Children');
children = children([1:end-3]);

handles = [UserData.plots.handle];
for i=1:length(children)
    ind = find(handles==children(i));
    if isempty(ind)
        ind = length(UserData.plots)+1;
        UserData.plots(ind).handle = children(i);
        UserData.plots(ind).pos0 = get(children(i),'Position');
    end
end

for i=1:length(UserData.plots);
    h = UserData.plots(i).handle;
    pos0 = UserData.plots(i).pos0;

    if length(pos0)<4
        continue;
    end
    pos(1) = -valueH+pos0(1)*valueZoom;
    pos(2) = 1+valueV - (valueZoom) + pos0(2)*valueZoom;

    if length(pos0)>2
        pos(3) = pos0(3) * valueZoom;
        pos(4) = pos0(4) * valueZoom;
    end
    set(h,'Position',pos);
end

set(fig,'UserData',UserData');

%{
function slider_callbacks(hObject, eventdata,handles);
% fig = gcf;
fig = get(hObject,'Parent');
UserData = get(fig,'UserData');

valueH = get(UserData.sliderH,'Value');
valueV = get(UserData.sliderV,'Value');
valueZoom = get(UserData.sliderZoom,'Value');
valueZoom = (1+(1-valueZoom)*4);

children = get(gcf,'Children');
% children = children([1:end-3,end]);
children = children([1:end-3]);

for i=1:length(children)
UserData.plots(i).handle = children(i);
%     pos0 = get(children(i),'OuterPosition');
pos0 = get(children(i),'Position');
if ~isfield(UserData.plots,'pos0') || isempty(UserData.plots(i).pos0)
UserData.plots(i).pos0 = pos0;
end
%         isempty(UserData.plots(i).pos0
%     UserData.plots(i).pos0 = children(i);
end

% [valueH,valueV,valueZoom]
for i=1:length(UserData.plots);
h = UserData.plots(i).handle;
pos = UserData.plots(i).pos0;

if length(pos)<4
continue;
end
pos0(1) = -valueH+pos(1)*valueZoom;
pos0(2) = 2-valueV - (valueZoom) + pos(2)*valueZoom;
pos0(3) = pos(3) * valueZoom;
pos0(4) = pos(4) * valueZoom;

%     set(h,'OuterPosition',pos);
set(h,'Position',pos0);
end

set(fig,'UserData',UserData');
%}