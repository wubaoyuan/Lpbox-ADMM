function shAsgImgs(Fs, Pts, Egs, asgs, Ax)
% Show image assignment result.
% 
% Input
%   Fs      -  images, 1 x 2 (cell)
%   Pts     -  points, 1 x 2 (cell)
%   Egs     -  edges, 1 x 2 (cell)
%   asgs    -  assignments, 1 x mAlg (cell)
%   Ax      -  axis, rows x cols (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-19-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-15-2012

% dimension
mAlg = length(asgs);

parCor = st('cor', 'ln', 'mkSiz', 7, 'cls', {'y', 'b'});

box0s = imgBox(Fs);
[boxs, boxG, Lns, Sca] = boxCat('horz', box0s, 'gap', 0);
Pt2s = gphBoxNew(Pts, box0s, boxs, Sca);

asgT = asgs{1};

for iAlg = 1 : mAlg
    asg = asgs{iAlg};
    ax = Ax{iAlg};
    
    shImgBox(Fs, boxs, 'ax', ax);
    shGph(Pt2s, Egs);
    shCor(Pt2s, asg.X, asgT.X, parCor);
    shLnCor(Lns, boxG);
    axis off;
    
    if isfield(asg, 'acc')
        title(sprintf('%s: %.2f', asg.alg, asg.acc));
    else
        title(sprintf('%s', asg.alg));        
    end
end
