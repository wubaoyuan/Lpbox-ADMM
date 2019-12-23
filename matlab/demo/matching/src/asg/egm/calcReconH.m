function [hs, errs] = calcReconH(Pts, idx, c, eta)
% Calculates the reconstruction weights for each point given a point-set.
% 
% Input
%   Pts     -  point set, 1 x m (cell), d x k
%   idx     -  neighbor index, 1 x nVar
%   c       -  row id
%   eta   -  eta for reguralization, {1000}
%
% Output    
%   hs      -  weight, 1 x k
%   errs    -  errors, 1 x m
%           
% History   
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-05-2012

% dimension
m = length(Pts);
[d, k] = size(Pts{1});
prIn('calcReconH', 'm %d, k %d, d %d', m, k, d);

% calculate the convex combination
nCon = 1 + d * m;
nVar = length(idx);

% set up sum(w) = 1
L = ones(nCon, nVar);
b = ones(nCon, 1);

Pt0s = Pts;
%Pts = falNormSca(Pts, c);

% set up d rows to be neighbors
for i = 1 : m
    pHd = 1 + (i - 1) * d;
    pEd = pHd + d - 1;
    L(pHd : pEd, :) = Pts{i}(:, idx);
    b(pHd : pEd) = Pts{i}(:, c);
end

% solve the linear system: A * w = b
if m == 1
    [w, flag] = lsqr(L, b, 1e-12, 10);
else
%    w = optL1('matlab3', L, b, 'th', 0);
     w = optL1('matlab2', L, b, 'th', eta);     
%    w = optL1('cvx', L, b);
end

% visualize
isDeb = 0;
if isDeb
    [Pt1s, PtMe, PtCs] = falNormSca(Pt0s, c);
    Pt2 = [];
    Pt3 = [];
    for i = 1 : m
        Pti = Pt1s{i};
        ptc = Pti(:, c);
        Pti = Pti - repmat(ptc, 1, 29);
        Pt2 = [Pt2, Pti];
        
        Pti = Pt0s{i};
        ptc = Pti(:, c);
        Pti = Pti - repmat(ptc, 1, 29);
        Pt3 = [Pt3, Pti];
    end
    rows = 1; cols = 3;
    Ax = iniAx(12, rows, cols, [250 * rows, 250 * cols] * 2);
    set(gcf, 'currentaxes', Ax{1});
    plot(Pt3(1, :), Pt3(2, :), 'ro');
    axis ij;
    axis equal;
    
    set(gcf, 'currentaxes', Ax{2});
    plot(Pt2(1, :), Pt2(2, :), 'ro');    
    axis ij;
    axis equal;
    
    set(gcf, 'currentaxes', Ax{3});
    plot(PtMe(1, :), PtMe(2, :), 'ro');
    hold on;
    plot(Pt1s{1}(1, :), Pt1s{2}(2, :), 'bs');
    plot(Pt0s{1}(1, :), Pt0s{2}(2, :), 'g^');    
    plot(PtCs(1, :, 1), PtCs(2, :, 1), 'k+');
    axis ij;
    axis equal;
end

% store
hs = zeros(1, k);
hs(idx) = -w;
hs(c) = 1;

% error
errs = zeros(1, m);
for i = 1 : m
    errs(:, i) = sum(abs(hs * Pts{i}'));
end

prOut;
