function shAsgShps2(Pts, Egs, asgs, Ax)
% Show graphs as well as the correspondence matrix.
% 
% Input
%   Pts     -  
%   Egs
%   asgs    -  assignment, 1 x mAlg (cell)
%   Ax      -  axis, 1 x (mAlg + 1) (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-19-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-14-2012

% dimension
mAlg = length(asgs);

parCor = st('cor', 'col', 'mkSiz', 7, 'cls', {'y', 'b'});

for c = 1 : mAlg
    X = asgs{c}.X;
    [n1, n2] = size(X);
    idx1 = 1 : n1;
    [~, idx2] = asgX2Idx(X);

    if c == 1
        set(gcf, 'CurrentAxes', Ax{1});
        shGph(Pts(1), Egs(1), 'parMks', {st('mkSiz', 0)});
        shPtCol(Pts{1}, idx1, n1, parCor);
        %title('template');
        axis off;
    end

    set(gcf, 'CurrentAxes', Ax{c + 1});
    shGph(Pts(2), Egs(2), 'parMks', {st('mkSiz', 0)});
    shPtCol(Pts{2}, idx2, n1, parCor);
    %title(asgs{c}.alg);
    axis off;
end
