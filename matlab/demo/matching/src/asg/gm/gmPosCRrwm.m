function X = gmPosCRrwm(K, X0, par)
% Reweighted Random Walk Matching.
%
% Reference
%   M. Cho, J. Lee, and K. M. Lee. "Reweighted random walks for
%   graph matching", In ECCV, 2010.
%
% Input
%   K       -  affinity matrix, nn x nn (sparse)
%   X0      -  initial assignment, n1 x n2
%   par     -  parameter
%     c     -  {0.2}
%                1.0 : for no reweighted jump
%                0.2 : used in the original paper
%
% Output
%   X       -  permutation matrix, n1 x n2
%
% History
%   create  -  Minsu Cho (chominsu@gmail.com), 09-25-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 11-03-2012

% function parameter
c = ps(par, 'c', .2);
isDisplay = 0;             % flag for visualization
amp_max = 30;              % maximum value for amplification procedure
iterMax = 300;             % maximum iterations for random walks 
thresConvergence = 1e-25;  % convergence threshold for random walks
tolC = 1e-3;               % convergence threshold for the Sinkhorn method

% dimension
[n1, n2] = size(X0);
E12 = ones(n1, n2);
nSize = n1 * n2;

% load data
M = K;

% transfrom E12 to list
[L12(:, 2), L12(:, 1)] = find(E12');

% eliminate conflicting elements to prevent conflicting walks
for i = 1 : nSize
    for j = i + 1 : nSize
        if L12(i, 1) == L12(j, 1) || L12(i, 2) == L12(j, 2) % one-to-one
            M(i, j) = 0; 
            M(j, i) = 0;
        end
    end
end

% make a transition matrix for affinity-preserving RW
options.disp = 0;
M = sparse(M);

% note that this matrix is column-wise stochastic
d = sum(M, 1); % degree : col sum
maxD = max(d);
Mo = M ./ maxD; % nomalize by the max degree

% initialize answer
prev_score = X0(:) / nSize;           % buffer for the previous score
prev_score2 = prev_score;             % buffer for the two iteration ahead
prev_assign = ones(nSize, 1) / nSize; % buffer for Sinkhorn result of prev_score

bCont = 1;
iter_i = 0;

if isDisplay
    figure(1);
end

% start main iteration
while bCont && iter_i < iterMax

    iter_i = iter_i + 1;

    if isDisplay 
        fprintf('iter %d\n', iter_i);
        subplot(231);
        imagesc(reshape(prev_score, n1, n2));
        title(sprintf('(%d) prev score',iter_i));
    end
    
    % random walking with reweighted jumps
    cur_score = Mo * (c * prev_score + (1 - c) * prev_assign);
    
    % normalization of sum 1
    sumCurScore = sum(cur_score);
    if sumCurScore > 0
        cur_score = cur_score ./ sumCurScore; 
    end
    
    if isDisplay
        subplot(232);
        imagesc(reshape(cur_score, n1, n2));
        title(sprintf('(%d) updated score', iter_i));
    end

    % update reweighted jumps
    cur_assign = cur_score;
    
    % attenuate small values and amplify large values
    amp_value = amp_max / max(cur_assign);     % compute amplification factor
    cur_assign = exp(amp_value * cur_assign);  % amplification 

    if isDisplay
        subplot(233);
        imagesc(reshape(cur_assign, n1, n2));
        title(sprintf('(%d) exponentiation',iter_i));
    end

    % Sinkhorn method of iterative bistocastic normalizations
    X = bistocNormalize_slack(reshape(cur_assign, n1, n2), tolC);
    cur_assign = X(:);
    
    sumCurAssign = sum(cur_assign); % normalization of sum 1
    if sumCurAssign > 0
        cur_assign = cur_assign ./ sumCurAssign; 
    end
    
    if isDisplay
        subplot(234);   
        imagesc(reshape(cur_assign, n1, n2));
        %bar3(reshape(cur_score, n1, n2),'detachted');
        title(sprintf('(%d) softassign',iter_i));
        %seeds = zeros(nSize,1);    seeds(find(answer)) = cur_assign(find(answer));
        subplot(236);
%        imagesc(problem.T12);
        %bar3(problem.T12,'detachted');
        title(sprintf('ground truth',iter_i));
        drawnow;
        pause(1);
        %pause;
    end
    
    % Check the convergence of random walks
    diff1 = sum((cur_score - prev_score) .^ 2);
    diff2 = sum((cur_score - prev_score2) .^ 2); % to prevent oscillations
    diff_min = min(diff1, diff2);
    if diff_min < thresConvergence
        bCont = 0;
    end

    prev_score2 = prev_score;
    prev_score = cur_score;
    prev_assign = cur_assign;
end

score = cur_score;
X = reshape(cur_assign, size(E12, 1), size(E12, 2));
