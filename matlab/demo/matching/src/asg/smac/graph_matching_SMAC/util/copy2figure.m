function copy2figure
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

currentFigure = gcf;
currentAxes = gca;
% gca
% posXLim = get(gca,'XLim');
% posYLim = get(gca,'YLim');
% q = posXLim(2)-posXLim(1);
% p = posYLim(2)-posYLim(1);
% % size(getimage);

if ~isempty(getimage(imgca))
    [p,q,r] = size(getimage(imgca));
% if ~isempty(getimage(gca))
%     [p,q,r] = size(getimage(gca));
else
    posXLim = get(gca,'XLim');
    posYLim = get(gca,'YLim');
    q = posXLim(2)-posXLim(1);
    p = posYLim(2)-posYLim(1);
end

figure2(get(currentFigure,'Name'),[p,q]);
newFigure = gcf;

newAxes = copyobj(currentAxes,newFigure);

figure(newFigure);
subplot(newAxes);

set(newAxes,'Position',[0 0 1 1]);
set(newFigure,'Colormap',get(currentFigure,'Colormap'));
% positionFigure;
% axis image;
% axis normal;
axis off;
% axis image;

positionFigure;%(p,q,newFigure);

% axis fill;
% gca
% newAxes
% z=get(gcf,'Children')
% positionFigure;
%     temp = get(currentAxes,'Children');
%     image = get(temp(length(temp)),'CData');

%     options = get(currentAxes,'UserData');
%     if ~isempty(options) && options.grid
%         afficheGrid(options.imageX,0,options.dataPatches.tableau_i,options.dataPatches.tableau_j);
%     end

