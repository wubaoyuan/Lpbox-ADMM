function [anchors, obj] = laicFindAnc(X, H, idxNeis, coeNeis, T, C, Reg)
% This function performs the "consistent-rounding" process described in
% Section II.B of H. Jiang, M. S. Drew, and Z. Li, "Matching by Linear 
% Programming and Successive Convexification", TPAMI 2007.
%
% Input
%   X        -  The Nm x Nt matching matrix. it does not need to be a binary matrix
%   H           The Nm x Nm reconstruction matrix of the template point set.
%   idxNeis  -  An N x 1 cell vector. The ith element records the ith point's 
%               neighbors' indices.
%   coeNeis  -  An N x 1 cell vector. The ith element records the ith point's 
%               neighbors' reconstruction weights.
%   T        -  An Nt x 2 matrix recording 2D Nt target points' coordinates.
%   C        -  The Nm x Nt feature matching cost matrix.
%   Reg      -  trust region, Nm x 4
%
% Output        
%   anchors  -  The Nm x 1 anchor vector. The ith element represents the
%               ith model point's corresponding target point index.
%   obj      -  Objective function's value anchors' result.
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 12-16-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 07-28-2012

prIn('laic find anchor');

% dimension
[Nm, Nt] = size(X);

% find related points
relatedT = cell(Nm, 1);
for i = 1 : Nm
    relatedT{i} = find(H(:, i) ~= 1 & H(:, i) ~= 0);
end
unrelatedT = H == 0;

% calculate each mapped feature cost
XFeatureCost = sum(C .* full(X), 2);

% calculate mapped points
MT = X * T;

% calculated mapped points geometric costs
% each row store one point's geometric cost
XgeoCost = sum(abs(H * MT), 2);

obj = 1e+30;
anchors = -1 * ones(Nm, 1);
% for each point in M find its anchor
for i = 1 : Nm
    minEn = 1e+30;
    
    % calculate feature cost subtotal for accelration
    ind = true(Nm, 1);
    ind(i) = false;
    XFeatureSubSum = sum(XFeatureCost(ind));
    
    % calculate geometric cost subtotal that is not related to current point
    XGeoSubSum = sum(XgeoCost(unrelatedT(:, i)));

    for j = 1 : Nt
        x = T(j, 1);
        y = T(j, 2);
        
        if x < Reg(i, 1) || x > Reg(i, 2) || y < Reg(i, 3) || y > Reg(i, 4)
            continue;
        end
        
        % feature cost
        featureCost = C(i, j);

        % leave one out, replace it with the tested one
        featureCost = featureCost + XFeatureSubSum;

        % geometric cost
        % current point geometric cost
        geoCost = sum(abs(T(j, :) - sum(MT(idxNeis{i}, :) .* repmat(coeNeis{i}, 1, 2))));

        % related points' geometric cost
        % this point might be related to several point's reconstruction
        for k = 1 : length(relatedT{i})
            relatedInd = relatedT{i}(k);

            selfInd = idxNeis{relatedInd} == i;
            otherInd = ~selfInd;
            
            geoCost = geoCost + sum(abs(MT(relatedInd, :) - T(j, :) * coeNeis{relatedInd}(selfInd) ...
                - sum(MT(idxNeis{relatedInd}(otherInd), :) .* repmat(coeNeis{relatedInd}(otherInd), 1, 2))));
        end

        % all other cost
        geoCost = geoCost + XGeoSubSum;
        
        % objective function value
        currentEn = featureCost + geoCost;
        if currentEn < minEn
            minEn = currentEn;
            anchors(i) = j;
            if minEn < obj
                obj = minEn;
            end
        end
    end
end

prOut;