function [PtD, dsts, angSs, angAs] = gphEg2Feat(Pt, Eg)
% Compute graph edge feature.
%
% Input
%   Pt      -  graph node, d x n
%   Eg      -  graph edge, 2 x m | []
%
% Output
%   PtD     -  edge vector, d x m
%   dsts    -  distance, 1 x m
%   angSs   -  angle (symmetric), 1 x m
%   angAs   -  angle (asymmetric), 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-26-2013

if isempty(Eg)
    PtD = [];
    dsts = [];
    angs = [];
    angAs = [];
    return;
end

% edge vector
Pt1 = Pt(:, Eg(1, :));
Pt2 = Pt(:, Eg(2, :));
PtD = Pt1 - Pt2;

% distance
dsts = real(sqrt(sum(PtD .^ 2)));

% angle (symmetric)
angSs = atan(PtD(2, :) ./ (PtD(1, :) + eps));

% angle (asymmetric)
angAs = atan2(PtD(2, :), PtD(1, :));
