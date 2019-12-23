function [X, best_score] = gmPosDIpfp(K, Ct, X0, par)
% This function tries to maximize the matching score x' M x 
% where x obeys discrete one-to-one matching constraints 
% such that x(i) = 1 if nodes(i) is matched to labels(i) and 0 otherwise.
%
% Reference
%   M. Leordeanu and M. Hebert and R. Sukthankar, "An Integer Projected
%   Fixed Point Method for Graph Matching and MAP Inference", in NIPS, 2009
%
% Remark
%   nn = n1 x n2
%
% Input
%   K        -  affinity matrix, nn x nn (sparse)
%   Ct       -  constraints, n1 x n2
%   X0       -  initial assignment, n1 x n2
%   par      -  parameter
%     nItMa  -  #maximum iteration steps, {50}
%
% Output
%   X        -  permutation matrix, n1 x n2
%
% History
%   create   -  Marius Leordeanu (leordeanu@gmail.com), 02-25-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 07-14-2012

% function parameter
nItMa = ps(par, 'nItMa', 50);
isDeb = psY(par, 'deb', 'n');

% dimension
[n1, n2] = size(X0);

M = K;
v = X0(:);

nNodes = n1;
nLabels = n2;
el = 0;
labels = zeros(1, nNodes * nLabels);
nodes = zeros(1, nNodes * nLabels);
for node = 1 : nNodes
    for label = 1 : nLabels
       el = el + 1;
       nodes(el) = node;
       labels(el) = label;
    end
end

% dimension
n = length(labels);
nLabels = max(labels);
nNodes = max(nodes);

% double-stochastic normalization
% v0 = v;
% v1 = v;
% for k = 1 : 20
%     for j = 1 : nNodes
%         f = find(nodes == j);
%         v1(f) = v0(f) / (sum(v0(f)) + eps);
%     end
%     for j = 1 : nLabels
%         f = find(labels == j);
%         v0(f) = v1(f) / (sum(v1(f)) + eps);
%     end
% end
% v = (v1 + v0) / 2;

if isDeb
    rows = 2; cols = 4;
    Ax = iniAx(9, rows, cols, [250 * rows, 250 * cols], 'pos', [0 0 .8 1]);
    ha = [];
end

% initial
sol0 = v;
new_sol = sol0;
best_sol = sol0;
best_score = 0; 
nSteps = 0;
scores(1) = new_sol' * M * new_sol;
scores2(1) = new_sol' * M * new_sol;
discreteRate = 0;
while nSteps <= nItMa
   nSteps = nSteps + 1;

   old_sol = new_sol;
   xx = M * old_sol;
   A = -Inf * ones(nNodes, nLabels);
   for i = 1 : nNodes
       f = find(nodes == i);
       A(i, labels(f)) = xx(f);
   end

   % gradient direction
   A = max(A(:)) - A;
   [X, score] = hungarian(A);
   x2 = zeros(n, 1);
   for i = 1 : nNodes
      f = find(nodes == i);
      match_ind = find(X(i, labels(f)) == 1);
      x2(f(match_ind)) = 1;
   end

   % step size
   k = (x2 - old_sol)' * M * (x2 - old_sol);
   t = 1;
   if k >= 0
       new_sol = x2;
       stepSize_t(nSteps) = 1;
       stepSize_norm(nSteps) = norm(x2 - old_sol);
       discreteRate = discreteRate + 1;
   else
       c = old_sol' * M * (x2 - old_sol);
       t = min([1, -c / k]);
       if t < 0.01
           t = 0;
       elseif t == 1
          discreteRate = discreteRate + 1;
       end
       new_sol = old_sol + t * (x2 - old_sol);
       stepSize_t(nSteps) = t;
       stepSize_norm(nSteps) = norm(x2 - old_sol);
   end

   scores(nSteps + 1) = new_sol' * M * new_sol;
   scores2(nSteps + 1) = new_sol' * M * old_sol;
   dX(nSteps) = sum(abs(new_sol - best_sol));
   curr_score = x2' * M * x2;
   
   if curr_score > best_score
        best_score = curr_score;
        best_sol = x2;
   end
   
   if isDeb
       ha = deb(ha, Ax, nSteps, scores, scores2, stepSize_t, best_sol, new_sol, x2, [n1 n2]);
   end
   
   % stop condition
   if norm(new_sol - old_sol) == 0
       break;
   end
end

discreteRate = discreteRate / nSteps;
sol = best_sol;
stats.dX = dX;
stats.scores = scores;
stats.scores2 = scores2;
stats.best_score = best_score;
stats.discreteRate = discreteRate;
stats.stepSize_t = stepSize_t;
stats.stepSize_norm = stepSize_norm;

% reshape
X = reshape(sol, [n1 n2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = deb(ha, Ax, nIt, scores, scores2, stepSize_t, best_sol, new_sol, x2, ns)
% Debug
%
% Input

% score
shIt(scores(1 : nIt), ones(1, nIt), 'ax', Ax{1, 1}, 'mkSiz', 7, 'itMa', 0); 
title('scores');
shIt(scores2(1 : nIt), ones(1, nIt), 'ax', Ax{2, 1}, 'mkSiz', 7, 'itMa', 0);
title('scores2');
shIt(stepSize_t(1 : nIt), ones(1, nIt), 'ax', Ax{1, 4}, 'mkSiz', 7, 'itMa', 0);
title('stepSize_t');

% new_sol
XNew = reshape(new_sol, ns);
if nIt == 1
    ha.hXNew = shM(XNew, 'ax', Ax{1, 2});
else
    shMUpd(ha.hXNew, XNew);
end

% x2
X2 = reshape(x2, ns);
if nIt == 1
    ha.hX2 = shM(X2, 'ax', Ax{1, 3});
else
    shMUpd(ha.hX2, X2);
end

% best_sol
XBest = reshape(best_sol, ns);
if nIt == 1
    ha.hXBest = shM(XBest, 'ax', Ax{2, 2});
else
    shMUpd(ha.hXBest, XBest);
end

if stepSize_t(nIt) > 1
    error('invalide step size');
end

drawnow;
pause(1);
