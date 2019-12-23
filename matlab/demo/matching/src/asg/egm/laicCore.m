function [X, C, H] = laicCore(Pts, E, C, lam)
% Core implementation of LAIC algorithm.
%
% Input
%   Pts     -  point set, 1 x 2 (cell), 2 x ni
%   E       -  graph adjacency matrix for source, Nm x Nm
%   C       -  feature matching cost matrix, Nm x Nt
%   lam     -  The parameter weighting feature cost and geometric cost, {1}
%
% Output
%   X       -  correspondence Hmatrix, Nm x Nt
%   H       -  LLE matrix, Nm x Nm
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-27-2012

% points
M = Pts{1}';
T = Pts{2}';

% dimension
[Nm, Nt] = size(C);

% apply lam to C
C = C ./ lam;

% calculate the reconstruction matrix
H = calcReconCoe(M, E);

% basis
bas = laicHull(T, C);

% construct LP0 and solve it
[X, obj] = laicLP(H, T, C, bas);
