function [sol, stats] = spectral_matching_ipfp(M, labels, nodes)
%
% Utility:
% this function tries to maximize the matching score x' M x 
% where x obeys discrete one-to-one matching constraints 
% such that x(i) = 1 if nodes(i) is matched to labels(i) and 0 otherwise
%
% History
%   create  -  Marius Leordeanu (leordeanu@gmail.com), 02-25-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-15-2011

% function parameter
maxSteps_IPFP = 50;
iterClimb_eigen = 30;

% dimension
n = length(labels);
nNodes = max(nodes);
nLabels = max(labels);

% compute the first eigenvector
v = ones(length(nodes), 1);
v = v / norm(v);
for i = 1 : iterClimb_eigen
  v = M * v;
  v = v / norm(v);
end

% double-stochastic normalization
v0 = v;
v1 = v;
for k = 1:20
    for j = 1 : nNodes
        f = find(nodes == j);
        v1(f) = v0(f) / (sum(v0(f)) + eps);
    end

    for j = 1 : nLabels
        f = find(labels == j);
        v0(f) = v1(f) / (sum(v1(f)) + eps);
    end
end
v = (v1 + v0) / 2;

% IPFP
sol0 = v;
new_sol = sol0;
best_sol = sol0;
best_score = 0;
nSteps = 0;
scores(1) = new_sol' * M * new_sol;
scores2(1) = new_sol' * M * new_sol;
discreteRate = 0;
converged = 0;
while 1
   if nSteps > maxSteps_IPFP
       break
   end

   nSteps =  nSteps + 1;
   old_sol = new_sol;
   xx = M * old_sol;

   A = -Inf * ones(nNodes, nLabels);
   for i = 1 : nNodes
       f = find(nodes == i);
       A(i, labels(f)) = xx(f);
   end
   
   A = max(A(:)) - A;
   [X, score] = hungarian(A);
   x2 = zeros(n, 1);
   
   for i = 1:nNodes
      f = find(nodes == i);
      match_ind = find(X(i, labels(f)) == 1);
      x2(f(match_ind)) = 1;
   end
       
   k = (x2 - old_sol)' * M * (x2 - old_sol);
   t = 1;
   if k >= 0
       new_sol = x2;
       stepSize_t(nSteps) = 1;
       stepSize_norm(nSteps) = norm(x2 - old_sol);
       if converged == 0
           discreteRate = discreteRate + 1;
       end
   else
       c = old_sol' * M * (x2 - old_sol);
       t = min([1, - c / k]);
       if t < 0.01
           t = 0;
       end
       if t == 1 && converged == 0
          discreteRate = discreteRate + 1;
       end
       
       new_sol = old_sol + t * (x2 - old_sol);
       stepSize_t(nSteps) = -c / k;
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
   if norm(new_sol - old_sol) == 0
       break;
   end
end

discreteRate = discreteRate / nSteps;

sol = best_sol;

% store
stats.dX = dX;
stats.scores = scores;
stats.scores2 = scores2;
stats.best_score = best_score;
stats.discreteRate = discreteRate;
stats.stepSize_t = stepSize_t;
stats.stepSize_norm = stepSize_norm;
