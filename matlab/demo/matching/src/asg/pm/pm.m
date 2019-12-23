function asg = pm(K, KQ, gphs, asgT)
% Optimal soft graph matching. Use exp((-d .* d) ./ sigma) as a similarity
% between edges, where d is the difference between the edges weight.
%
% Reference
%   R. Zass and A. Shashua, "Probabilistic Graph and Hypergraph Matching",
%   in CVPR, 2008
%
% Input
%   K       -   global kernel matrix, nn x nn
%   KQ      -   edge kernel matrix, m1 x m2
%   gphs    -  graphs, 1 x 2 (cell)
%   asgT    -  ground-truth assignment (can be [])
%
% Output
%   asg
%
% History
%   create  -  Ron Zass (zass@cs.huji.ac.il), 01-25-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-26-2012

% function parameter
prIn('pm');

% graph element
G1 = gphs{1}.G;
G2 = gphs{2}.G;

% dimension
n = size(G1, 1);

% marginalize
Y = G1 * KQ * G2';

% hypergraph matching over Y
Z = hypergraphMatching(Y, n);

% post-processing
X = gmPosDHun(Z);

% compare with ground-truth
acc = matchAsg(X, asgT);

% store
asg.alg = 'pm';
asg.X = X;
asg.acc = acc;
asg.obj = X(:)' * K * X(:);

prOut;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = hypergraphMatching(Y, numberOfMatches)
% Optimal soft hyergraph matching.
%
% Input
% Y         - Marginalization of the hyperedge-to-hyperedge correspondences matrix.
% numberOfMatches - number of matches required.
%
% Output
%   Z       - an n1 by n2 matrix with the soft matching results.

% dimension
[n1, n2] = size(Y);

maxRowSum = ones(n1, 1);
maxColumnSum = ones(1, n2);

Z = nearestDSmax_RE(Y, maxRowSum, maxColumnSum, numberOfMatches);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = nearestDSmax_RE(Y, maxRowSum, maxColSum, totalSum, precision, maxLoops)
% Find the nearest matrix with non-negative entries and row / column sum at
% most rowSum / colSum, to A in relative entropy.
%
% Input
%   Y          - an m by n non-negative matrix.
%   maxRowSum  - a non-negative column vector of length m with the maximal allowed rows sum.
%   maxColSum  - a non-negative row vector of length n with the maximal allowed columns sum.
%   totalSum   - a non-negative scalar, the desired sum of all elements.
%   precision [Optional, default=0.01] - break when the each row/column sum is at most precision away from the desired sum.
%   maxLoops  [Optional, defalut=1000] - Max number of iterations to perform.
%
% Output
%   F          -  The nearest generalized doubly stochastic F to Ker in relative entropy.
    
if (nargin < 6)
    maxLoops = 1000;
    if (nargin < 5)
        precision = 0.01;
    end
end

lambda1 = zeros(size(Y));
lambda2 = lambda1;
lambda3 = lambda1;
F1 = totalSum * (Y ./ sum(Y(:)));
F2 = F1;
F3 = F1;

for t = 1 : maxLoops
    
    % Max row sum:
    H1 = lambda1 - (Y ./ (F3 + eps));
    F1 = maxColSumP(Y', -H1', maxRowSum', precision)';
    F1 = zeroNan(F1);
    lambda1 = lambda1 - (Y ./ (F3 + eps)) + (Y ./ (F1 + eps));
    
    % Max column sum:
    H2 = lambda2 - (Y ./ (F1 + eps));
    F2 = maxColSumP (Y, -H2, maxColSum, precision);
    F2 = zeroNan(F2);
    lambda2 = lambda2 - (Y ./ (F1 + eps)) + (Y ./ (F2 + eps));

    % Total sum:
    H3 = lambda3 - (Y ./ (F2 + eps));
    F3 = reshape(exactTotalSum(Y(:), -H3(:), totalSum, precision), size(Y));
    F3 = zeroNan(F3);
    lambda3 = lambda3 - (Y ./ (F2 + eps)) + (Y ./ (F3 + eps));

    if max(abs(F1(:) - F2(:))) < precision && max(abs(F1(:) - F3(:))) < precision
        break;
    end
end

F = (F1 + F2 + F3) / 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = unconstrainedP(Y, H)
X = Y ./ H;
X(find(X < eps)) = eps;

%%%%%%%%%%%%%%%%%%%%%%%
function F = zeroNan(F)

idx = isnan(F);
F(idx) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = exactTotalSum(y, h, totalSum, precision)
% y and h are vectors, totalSum and precision are scalars

totalSumMinus = totalSum - precision;
curAlpha = -min(h) + eps;

stepAlpha = max(10, abs(curAlpha / 10));
for i = 1 : 50 % inner maxLoops
    newAlpha = curAlpha + stepAlpha;
    x = y ./ (h + newAlpha + eps);
    newSum = sum(x);
    if newSum > totalSum
        curAlpha = newAlpha;
    else
        if newSum < totalSumMinus
            stepAlpha = stepAlpha / 2;
        else % sum between totalSumMinus and totalSum
            return;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = maxColSumP(Y, H, maxColSum, precision)
X = unconstrainedP(Y, H);
Xsum = sum(X);
for i = find(Xsum > maxColSum)
    X(:, i) = exactTotalSum(Y(:,i), H(:,i), maxColSum(i), precision);
end
