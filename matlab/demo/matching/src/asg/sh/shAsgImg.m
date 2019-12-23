function shAsgImg(Fs, gphs, asg, asgT, parCor, varargin)
% Show image with correspondence.
% 
% Input
%   Fs      -  images, 1 x 2 (cell)
%   gphs    -  graphs, 1 x 2 (cell)
%   asg     -  assignment
%   parCor  -  parameter for showing correspondence
%   varargin
%     show option
%     ord   -  flag of re-order points, {'y'} | 'n'
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-19-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-02-2013

% show option
psSh(varargin);

% function option
isOrd = psY(varargin, 'ord', 'n');

% graph nodes
Pts = {gphs{1}.Pt, gphs{2}.Pt};
Egs = {gphs{1}.Eg, gphs{2}.Eg};

% graph plotting property
parMks = {[], []};

% re-order
if isOrd
    [Pts, Egs, asg] = gphReOrd(Pts, Egs, asg);
end

% concate images
box0s = imgBox(Fs);
[boxs, boxG, Lns, Sca] = boxCat('horz', box0s, 'gap', 0);
Pt2s = gphBoxNew(gphs, box0s, boxs, Sca);

% plot image
shImgBox(Fs, boxs);

% plot graph nodes and edges
shGph(Pt2s, Egs);

% plot correspondence
XT = ps(asgT, 'X', []);
shGphCor(Pt2s, asg.X, XT, parCor);

% plot line for splitting
shBreakLn(Lns, boxG);

axis off;
