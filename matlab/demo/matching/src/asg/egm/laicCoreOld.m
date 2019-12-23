function [X, Xs, C, H, idxTrs, idxTes] = laicCoreOld(Pts, E, C, par)
% Core implementation of LAIC algorithm.
%
% References
%   H. Li, E. Kim, X. Huang, and L. He, 
%   "Object Matching Using a Locally Affine-Invairant Constraint", in CVPR, 2010.
%
% Input
%   Pts       -  point set, 1 x 2 (cell), 2 x ni
%   E         -  graph adjacency matrix for source, Nm x Nm
%   C         -  feature matching cost matrix, Nm x Nt
%   par
%     lambda  -  The parameter weighting feature cost and geometric cost, {1}
%     nItR    -  The function will use the the "consitent-rounding" for the last nItR iterations, {2}
%                see Secion II.B of Jiang et al, TPAMI 2007.
%     dias    -  should be [maxDiameter minDiameter diameterSpeed finalDiameter], {[150 15 15 30]}
%                maxDiameter:   the maximum diameter of trust regions.
%                minDiameter:   the minimum diameter of trust regions.
%                diameterSpeed: the diameter decreasing speed.
%                finalDiameter: the final diameter for non-simplified LP model.
%     Nsame   -  A parameter tuning the "matching to the same target point"
%                constraint. For each target point, no more than Nsame model
%                points can be matched to it, {3}
%
% Output
%   X         -  correspondence Hmatrix, Nm x Nt
%   H         -  LLE matrix, Nm x Nm
%
% History
%   create    -  Hongsheng Li (h.li@lehigh.edu), 12-13-2010
%   modify    -  Feng Zhou (zhfe99@gmail.com), 07-27-2012

% function parameters
lambda = ps(par, 'lambda', 1);
nItR = ps(par, 'nItR', 2);
dias = ps(par, 'dias', [150 15 15 30]);
Nsame = ps(par, 'Nsame', 3);
prIn('laic', 'lambda %.0f, nItR %d, Nsam %d, dias [%d %d %d %d]', lambda, nItR, Nsame, dias(1), dias(2), dias(3), dias(4));

% points
M = Pts{1}';
T = Pts{2}';

% dimension
Nm = size(M, 1);
Nt = size(T, 1);

Xs = cell(1, 100);
mX = 0;

% get diameter parameters
maxDiameter = dias(1);
minDiameter = dias(2);
diameterSpeed = dias(3);
finalDiameter = dias(4);

% #iteration
iterNum = floor((maxDiameter - minDiameter) / diameterSpeed);
iterNumThres = iterNum - nItR;
if nItR > iterNum || nItR < 1
    error('nItR should be between 1 and iterNum(%d)!', iterNum);
end

% apply lambda to C
C = C ./ lambda;

% calculate the reconstruction matrix
[H, idxNeis, coeNeis] = calcReconCoe(M, E);

% initial basis
basis = cell(Nm, 1);
for i = 1 : Nm
    basis{i} = 1 : Nt;

    % convex hull
    costSurf = [T(basis{i}, :), C(i, basis{i})'];
    costSurfHull = convhulln(costSurf, {'Qt' 'PD2'});

    % support points of the convex hull
    costSurfHull = sort(costSurfHull(:));
    ind = costSurfHull ~= [costSurfHull(2 : end); 0];
    basis{i} = costSurfHull(ind);
    
    % visualize by feng
%     figure(10); clf;
%     trisurf(costSurfHull, costSurf(:, 1), costSurf(:, 2), costSurf(:, 3));
%     hold on;
%     plot3(costSurf(:, 1), costSurf(:, 2), costSurf(:, 3), 'ro', 'MarkerSize', 10);
%     plot3(costSurf(basis{1}, 1), costSurf(basis{1}, 2), costSurf(basis{1}, 3), 'b+', 'MarkerSize', 10);
end

% construct LP0 and solve it
[X, obj] = laicLPOld(H, T, C, basis, Nsame);
mX = mX + 1;
Xs{mX} = X;

% calculate anchor points using LP's result
anchorPoints = X * T;

% update trust regions
trustRegion = laicRegSet(anchorPoints, maxDiameter);

% shrink the trust region
for iter = 1 : 0
    pr('iter %d, diameter: %f', iter, maxDiameter - diameterSpeed * iter);
    beginTime = cputime;

    % for every point in the M, update trust region and update basis.
    for i = 1 : Nm
        
        % check whether each model point is inside its trust region
        insideInd = [];
        for j = 1 : Nt
            x = T(j, 1);
            y = T(j, 2);

            if x >= trustRegion(i, 1) && x <= trustRegion(i, 2) && y >= trustRegion(i, 3) && y <= trustRegion(i, 4)
                insideInd = [insideInd, j];
            end
        end

        % update basis by re-constructing the convex hull
        if length(insideInd) >= 4
            costSurf = [T(insideInd, :) C(i, insideInd)'];
            costSurfHull = convhulln(costSurf, {'Qt' 'PD2'});
            costSurfHull = sort(costSurfHull(:));
            ind = costSurfHull ~= [costSurfHull(2 : end); 0];
            basis{i} = insideInd(costSurfHull(ind));
        end
    end

    % construct LP and solve it
    [X, obj] = laicLPOld(H, T, C, basis, Nsame);
    mX = mX + 1;
    Xs{mX} = X;

    if iter <= iterNumThres
        % calculate anchor points directly using LP's result
        anchorPoints = X * T;

    else
        [anchors, newIntObj] = laicFindAnc(X, H, idxNeis, coeNeis, T, C, trustRegion);

        for ptInd = 1 : length(anchors)
            anchorPoints(ptInd, :) = T(anchors(ptInd), :);
        end
    end

    % update trust regions
    newTrustRegion = laicRegUpd(anchorPoints, trustRegion, diameterSpeed);

    emptyInd = laicRegEmp(newTrustRegion, T);
    trustRegion(~emptyInd, :) = newTrustRegion(~emptyInd, :);

    endTime = cputime;
    pr('iter %d, %f seconds', iter, endTime - beginTime);
end

% the last iteration without the convex trick
for i = 1 : Nm
     trustRegion = laicRegSet(anchorPoints, finalDiameter);

     insideInd = [];
     for j = 1 : Nt
         x = T(j, 1);
         y = T(j, 2);

         if x >= trustRegion(i, 1) && x <= trustRegion(i, 2) && y >= trustRegion(i, 3) && y <= trustRegion(i, 4)
             insideInd = [insideInd j];
         end
     end

     basis{i} = insideInd;
end

% resolve the un-simplified LP model
[X, obj] = laicLPOld(H, T, C, basis, Nsame);
mX = mX + 1;
Xs{mX} = X;

[anchors, newIntObj] = laicFindAnc(X, H, idxNeis, coeNeis, T, C, trustRegion);
X = zeros(size(X));
for i = 1 : Nm
    X(i, anchors(i)) = 1;
end
mX = mX + 1;
Xs{mX} = X;

Xs(mX + 1 : end) = [];

% matching with ground-truth assignment if possible
%[idxTrs, idxTes] = laicSortX(X, Pts, C, H);

prOut;
