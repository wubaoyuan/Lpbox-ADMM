function shAsgShps(Pts, Egs, asgs, Ax)
% Show the correspondence matrix.
% 
% Input
%   Pts     -  graph node, 1 x k (cell), 2 x ni
%   Egs     -  graph edge, 1 x k (cell), 2 x mi | []
%   asgs    -  assignment, 1 x mAlg (cell)
%   Ax      -  axis, rows x cols (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-19-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-19-2012

% dimension
mAlg = length(asgs);

% bounding box
box0s = gphBox(Pts);
[boxs, boxG, Lns, Sca] = boxCat('horz', box0s, 'gap', .2);
Pt2s = gphBoxNew(Pts, box0s, boxs, Sca);
parCor = st('cor', 'col');

for iAlg = 1 : mAlg
    asg = asgs{iAlg};
    ax = Ax{iAlg};
    
    shGph(Pt2s, Egs, 'ax', ax);
    shCor(Pt2s, asg.X, [], parCor);
    shBreakLn(Lns, boxG);
    
    if isfield(asg, 'acc')
        title(sprintf('%s: %.2f', asg.alg, asg.acc));
    elseif isfield(asg, 'alg')
        title(sprintf('%s', asg.alg));
    end
end
