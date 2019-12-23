function KEg = knlGphQ2E(KQ)
% Convert from sparse graph affinity to node and edge affinity.
%
% Remark
%   k0 = k10 x k20
%   In most case, k10 = k1 and k20 = k2.
%
% Input
%   KQ      -  2nd order affinity, k0 x k0 (sparse)
%
% Output
%   KEg     -  edge affinity, l1 x l2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-29-2012

% dimension
k0 = size(KQ, 1);

% default visPts
if isempty(visPts)
    k10 = k1;
    k20 = k2;
    
    visPts = cell(1, 2);
    visPts{1} = ones(k1, 1);
    visPts{2} = ones(k2, 1);
    
else
    k10 = size(visPts{1}, 1);
    k20 = size(visPts{2}, 1);
end
k0 = k10 * k20;

% index
idxPt1 = find(visPts{1})';
idxPt2 = find(visPts{2})';

% convert node affinity
KP = zeros(k10, k20);
KP(idxPt1, idxPt2) = KPt;

% convert edge affinity
I11 = repmat(idxPt1(Egs{1}(1, :))', 1, l2);
I12 = repmat(idxPt1(Egs{1}(2, :))', 1, l2);
I21 = repmat(idxPt2(Egs{2}(1, :)), l1, 1);
I22 = repmat(idxPt2(Egs{2}(2, :)), l1, 1);
I1 = sub2ind([k1 k2], I11(:), I21(:));
I2 = sub2ind([k1 k2], I12(:), I22(:));
KQ = sparse(I1(:), I2(:), KEg(:), k0, k0);
