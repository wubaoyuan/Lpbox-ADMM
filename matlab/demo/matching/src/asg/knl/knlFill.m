function [KP, Egs] = knlFill(KP0, Eg0s, visPts)
% Filling the missing nodes.
%
% Input
%   KP0     -  original node affinity, n1 x n2
%   Eg0s    -  original graph edges, 1 x 2 (cell), 2 x mi
%   visPts  -  node existence status, [] | 1 x 2 (cell), niT x 1
%
% Output
%   KP      -  new node affinity, n1T x n2T
%   Egs     -  new graph edges, 1 x 2 (cell), 2 x mi
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
[n1, n2] = size(KP0);

% index of point
if isempty(visPts)
    n1T = n1;
    n2T = n2;
    idxPt1 = 1 : n1;
    idxPt2 = 1 : n2;

else
    n1T = size(visPts{1}, 1);
    n2T = size(visPts{2}, 1);
    idxPt1 = find(visPts{1})';
    idxPt2 = find(visPts{2})';
end

% convert node affinity
KP = zeros(n1T, n2T);
KP(idxPt1, idxPt2) = KP0;

% edge
Eg1 = [idxPt1(Eg0s{1}(1, :)); idxPt1(Eg0s{1}(2, :))];
Eg2 = [idxPt2(Eg0s{2}(1, :)); idxPt2(Eg0s{2}(2, :))];
Egs = {Eg1, Eg2};
