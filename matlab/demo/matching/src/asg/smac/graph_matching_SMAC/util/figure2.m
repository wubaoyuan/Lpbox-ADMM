function fig = figure2(name, position, options, fig)
% Timothee Cour, 21-Apr-2008 17:31:23

h=get(0,'CurrentFigure');
if ishandle(h)
    pos=get(gcf,'OuterPosition');
else
    pos=[];
end

if nargin < 4
    global data;
    if ~isfield(data,'gui')
        data=initializeData(data);
    end
    if ~isfield(data.gui,'newWindow')
        data=initializeData(data);
    end

    if data.gui.newWindow ==1
        while ishandle(data.gui.indexNextFigure)
            data.gui.indexNextFigure=data.gui.indexNextFigure+1;
        end
    else
        if data.gui.indexNextFigure < 1
            data.gui.indexNextFigure = 1;
        end
    end
    data.figureCourante = figure(data.gui.indexNextFigure);
    clf;
else
    figure(fig);
end


if nargin >= 1
    set(gcf,'Name',name);
end
try %in case no display available
    if nargin >=2
        positionFigure(position);
    else
        if isempty(pos)
            positionFigure(1,1);
        else
            set(gcf,'OuterPosition',pos);
        end
    end
%     [p0i1,p0i2,p0j1,p0j2] = position2lim(get(0,'ScreenSize'));
%     p0j1=p0j1+274;
%     setFigureLocation(p0i1,p0j1,gcf);
catch
end
if nargin >= 3
    set(gca,'UserData',options);
end
axis off

setPropertyFigure2(gcf);
% set(gcf,'WindowButtonDownFcn',@actionFigure);

fig = gcf;

global data;
data.gui.dataFigures(fig).isClosable=1;
