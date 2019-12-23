function shAsgShp(gphs, asg, XT, parMks, parCor, varargin)
% Show shape with correspondence.
% 
% Input
%   gphs    -  graphs, 1 x 2 (cell)
%   asg     -  assignment
%   parCor  -  parameter for showing correspondence
%   show option
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-19-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 11-20-2012

% show option
psSh(varargin);

% concate shapes
box0s = gphBox(gphs);
[boxs, boxG, Lns, Sca] = boxCat('horz', box0s, 'gap', .2);
gph2s = gphBoxNew(gphs, box0s, boxs, Sca);

hold on;
% plot correspondence
shGphCor(gph2s, asg.X, XT, parCor);

% plot graph nodes and edges
shGph(gph2s, 'parMks', parMks);

% plot line for splitting
shBreakLn(Lns, boxG);

% axis off;
