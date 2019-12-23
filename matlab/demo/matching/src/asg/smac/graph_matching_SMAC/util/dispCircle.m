function h=dispCircle(xy,R,color,width);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if nargin<3
    color='b';
end
if nargin<4
    width=1;
end

t=0:0.01:2*pi;
hold on;
% plot(R*cos(t)+xy(2),R*sin(t)+xy(1));
h=plot(R*cos(t)+xy(2),R*sin(t)+xy(1),'color',color,'LineWidth',width);
hold off
% axis on
