function positionFigure(p,q,fig);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

set(gca,'Position',[0 0 1 1]);
if nargin < 3
    fig = gcf;
end

if nargin < 1
    if ~isempty(getimage(imgca))
        [p,q,r] = size(getimage(imgca));
%     if ~isempty(getimage(gca))
%         [p,q,r] = size(getimage(gca));
    else
        posXLim = get(gca,'XLim');
        posYLim = get(gca,'YLim');
        q = posXLim(2)-posXLim(1);
        p = posYLim(2)-posYLim(1);
    end

%     [p,q,r]=size(getimage(gca));
%     posXLim = get(gca,'XLim');
%     posYLim = get(gca,'YLim');
    
%     q = posXLim(2)-posXLim(1);
%     p = posYLim(2)-posYLim(1);

    %     pq = size(getimage)
    %     %     pq = size(get(get(gca,'Children'),'cdata'));
    %     p = pq(1);
    %     q = pq(2);
elseif nargin < 2 && length(p) == 2
    q = p(2);
    p = p(1);
elseif nargin < 2 && length(p) == 1
    q = 1;
elseif nargin < 2 && length(p) == 4
    set(fig,'Position',p);
    q = p(2);
    p = p(1);
elseif nargin >=2 && strcmp(p,'width')
    widthAxis = q;
    posXLim = get(gca,'XLim');
    posYLim = get(gca,'YLim');
    q = posXLim(2)-posXLim(1);
    p = posYLim(2)-posYLim(1);
end

positionFigure2(p,q);

% widthAxis = 300;%300
% heightBottomBorder = 150;
% widthGui = 278;
% heightAxis = widthAxis * p / q;
% screensize = get(0,'ScreenSize');
% if screensize(4)-heightAxis < heightBottomBorder
%     heightAxis = screensize(4) - heightBottomBorder;
%     widthAxis = heightAxis/(p / q);
% end
% set(gca,'Position',[0 0 1 1]);
% % set(fig,'Position',[4, 100, widthAxis, heightAxis ]);
% set(fig,'Position',[widthGui, screensize(4)-heightAxis, widthAxis, heightAxis ]);
% [widthGui, screensize(4)-heightAxis, widthAxis, heightAxis ]
% posFig = get(gcf,'OuterPosition');
% posFig
% heightBorderWindow = 26;
% if posFig(2)+posFig(4) - screensize(4) > heightBorderWindow
%     posFig(2) = posFig(2) - (posFig(2)+posFig(4) - screensize(4));
%     set(fig,'OuterPosition',posFig);
% %     set(fig,'OuterPosition',[274   670   308   380]);
%     posFig
% end
% 
% axis off;


function positionFigure2(p,q);
[p0i1,p0i2,p0j1,p0j2] = position2lim(get(0,'ScreenSize'));

% [p,q,r]=size(getimage(gca));

widthAxis = 500;%300
% q = round(q*widthAxis/(p+eps));
% p = widthAxis;
ratio = p/q;
q = widthAxis;
p = round(widthAxis*ratio);

borderBottom = 200;
if p > p0i2-borderBottom
    p = p0i2-borderBottom;
    q = round(p/ratio);
end
% q = round(q*widthAxis/(p+eps));


[p1i1,p1i2,p1j1,p1j2] = position2lim(get(gcf,'Position'));
[p2i1,p2i2,p2j1,p2j2] = position2lim(get(gcf,'OuterPosition'));

heightTitleBar = (p2i2-p2i1)-(p1i2-p1i1);
% heightTitleBar
if heightTitleBar < 80 %should be 80
    heightTitleBar = 80;
end
heightTitleBar = 35;%when no title bar

% p3i1 = p0i1;
% p3j1 = 274;
p3i1 = p2i1;
p3j1 = p2j1;
p3i2 = p3i1 + p + heightTitleBar;
p3j2 = p3j1 + q + (p2j2-p2j1)-(p1j2-p1j1);

set(gcf,'OuterPosition',lim2position([p3i1,p3i2,p3j1,p3j2]));
return;

% problem : while typing on command line or from figure call-back, we get
% different results for one of those : 
% get(gcf,'OuterPosition')
% get(gcf,'Position')
% A=[p0i1,p0i2,p0j1,p0j2;p1i1,p1i2,p1j1,p1j2;p2i1,p2i2,p2j1,p2j2];
% disp(round(A));