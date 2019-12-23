function [sol, x_opt, score, score_sol]  = ipfp(M, D, sol0, labels, nodes, maxIter)
% author: Marius Leordeanu
% last updated: Feb 25, 2011
%
% for questions contact the author at: leordeanu@gmail.com
%
% please cite the following paper:
%  
% An Integer Projected Fixed Point Method for Graph Matching and MAP
% Inference, NIPS 2009
% by Marius Leordeanu ,  Martial Hebert ,  Rahul Sukthankar
%
% Utility:
% this function tries to maximize the matching score x'Mx + Dx 
% (note that D containts the unary terms and M containts the pairwise terms)
% where x obeys discrete one-to-one matching constraints 
% such that x(i) = 1 if nodes(i) is matched to labels(i) and 0 otherwise
%
% Note: sol0 is the initial solution

% initializations
n = length(labels);
nLabels = max(labels);
nNodes  = max(nodes);
xc = sol0;
sol = sol0;
score_sol = 0;
score(1) =  xc' * M * xc + D' * xc;

% climbing ------------------------------------------------------------------
for nSteps = 1 : maxIter          
   x0 = xc;
   xx = M * x0 + D / 2;
   A = -Inf * ones(nNodes, nLabels);
   for i = 1 : nNodes
       f = find(nodes == i);
       A(i, labels(f)) = xx(f);
   end
   A = max(A(:)) - A;
   [X, s_aux] = hungarian(A);
   b = zeros(n,1);
   for i = 1 : nNodes
      f = find(nodes == i);
      match_ind = find(X(i, labels(f)) == 1);
      b(f(match_ind)) = 1;
   end

   dx = b - x0;
   A = dx' * M * dx;
   t = 1;
   
   if A < 0
       C = (x0' * M + D' / 2) * dx;
       t = min([1, -C / A]);
       if t < 0.01
           t = 0;
       end
   end
   
   xc = x0 + t * dx;
   
   score(nSteps + 1) = xc' * M * xc + D' * xc;
     
   score_b = b' * M * b + D' * b;
   
   if score_b >= score_sol
       
        score_sol = score_b;
        
        sol = b;
 
   end
      
   if norm(xc - x0) == 0
       break;
   end
   
end

x_opt = xc;


