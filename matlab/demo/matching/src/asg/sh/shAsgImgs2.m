function shAsgImgs2(Fs, gphs, asgs, Ax, varargin)
% Show image with correspondence.
% 
% Input
%   Fs      -  images, 1 x 2 (cell)
%   gphs    -  graphs, 1 x 2 (cell)
%   asg     -  assignment
%   Ax      -  axes
%   varargin
%     ord   -  flag of re-order points, {'y'} | 'n'
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-19-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-11-2012

% function option
isOrd = psY(varargin, 'ord', 'y');

% dimension
mAlg = length(asgs);

% graph nodes
Pts = {gphs{1}.Pt, gphs{2}.Pt};
Egs = {gphs{1}.Eg, gphs{2}.Eg};

% parameter for plotting the correspondence
parCor = st('cor', 'col', 'mkSiz', 7, 'cls', {'y', 'b'});

% re-order
if isOrd
    [Pts, Egs, asgs] = gphReOrd(Pts, Egs, asgs);
end

for c = 1 : mAlg
    X = asgs{c}.X;
    [n1, n2] = size(X);
    idx1 = 1 : n1;
    [~, idx2] = asgX2Idx(X);

    if c == 1
        shImg(Fs{1}, 'ax', Ax{1});
        shGph(gphs(1), 'parMks', {st('mkSiz', 0, 'lnWid', 1)});
        shPtCol(Pts{1}, idx1, n1, parCor);
        %title('template');
        axis off;
    end

    shImg(Fs{2}, 'ax', Ax{c + 1});
    shGph(gphs(2), 'parMks', {st('mkSiz', 0)});
    shPtCol(Pts{2}, idx2, n1, parCor);
    %title(asgs{c}.alg);
    axis off;
end
