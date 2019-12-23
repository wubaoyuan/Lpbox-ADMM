function Ps = asgMBB(C0s)
% Exact solution of MDA by using a Branch-and-Bound (BB) algorithm.
%
% Input
%   C0s     -  confusion matrix before adjusting, k x k x n x n
%
% Output
%   Ps      -  matching matrix, k x k x n x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-03-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% global matrix that will be frequently updated in the branch-and-bound algorithm
global PGs;
global CGs;
global visG;
global n;
global m;     % current level
global idxs;
global corr;
global isHug; % flog to use Hungarain algorithm to find the best match
global maxV;
global maxPs;
global bounds;

[k, k, n, n] = size(C0s);

% multi-dimensional assignment
PGs = zeros(k, k, n, n); CGs = C0s; visG = ones(k, n);
maxV = -inf; maxPs = []; m = k;

maxState = k ^ n;
bounds = zeros(m, maxState);
corr = zeros(m, n);
idxs = cell(m, n);

isHug = false;

tic;
treeSearch(0);
t = toc;
prom('m', 'Branch-and-Bround costs %.2f seconds for k = %d n = %d\n', t, k, n);
Ps = maxPs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = treeSearch(v0)
% Tree-like searching for a branch-and-bound algorithm.
%
% Input
%   v0  -  the objective value that is cumulatively calculated from the past path
%
% Output
%   v   -  maximum objective achieved in the current sub-tree (current level to the leaves)

global PGs;
global CGs;
global visG;
global n;
global m;
global idxs;
global corr;
global isHug;
global maxV;
global maxPs;
global bounds;

% bottom level, completed
if m == 0
    v = 0;
    if v0 > maxV
        maxV = v0;
        maxPs = PGs;
    end

    return;
end

% non-filled class
for i = 1 : n
    idxs{m, i} = find(visG(:, i));
end

% enumerate all possible correspondences
nState = m ^ (n - 1);
i = 1;
boundMiG = -inf;
while i <= nState

    corr(m, 1) = 1;
    i0 = i - 1;    
    for j = 2 : n
        corr(m, j) = rem(i0, m) + 1;
        i0 = floor(i0 / m + eps);
    end
    
    % obtain the upper bound of each possible correspondence
    boundMa = 0;
    boundMi = 0;
    for ii = 1 : n
        idxi = idxs{m, ii};
        ci = idxi(corr(m, ii));
        idxi(corr(m, ii)) = [];
        for jj = ii + 1 : n
            idxj = idxs{m, jj};
            cj = idxj(corr(m, jj));
            idxj(corr(m, jj)) = [];            

            boundMa = boundMa + CGs(ci, cj, ii, jj);
            boundMi = boundMi + CGs(ci, cj, ii, jj);

            C = CGs(idxi, idxj, ii, jj);

            if isHug
                if size(C, 1) > 1
                    P = assign(C);
                    val = trace(C * P');
                elseif size(C, 1) == 1
                    val = C;
                else
                    val = 0;
                end
                bound = bound + val;
            else
                valMa = max(C);
                valMi = min(C);
                
                boundMa = boundMa + sum(valMa);
                boundMi = boundMi + sum(valMi);                
            end
        end
    end
    bounds(m, i) = boundMa;
    
    if boundMi > boundMiG
        boundMiG = boundMi;
    end
    
    i = i + 1;
end

% sort the correspondence based on the maximum bound
% [boundsSort, ind] = sort(bounds(m, 1 : nState), 'descend');

% search each correspondence
minBound = boundMiG;
i = 1; 
while i <= nState

    if bounds(m, i) < minBound
        i = i + 1;
        continue;
    end

    corr(m, 1) = 1;
    i0 = i - 1;    
    for j = 2 : n
        corr(m, j) = rem(i0, m) + 1;
        i0 = floor(i0 / m + eps);
    end

    % fill
    vCurrent = 0;
    for ii = 1 : n
        ci = idxs{m, ii}(corr(m, ii));
        visG(ci, ii) = 0;
        
        for jj = ii + 1 : n
            cj = idxs{m, jj}(corr(m, jj));
            
            vCurrent = vCurrent + CGs(ci, cj, ii, jj);
            PGs(ci, cj, ii, jj) = 1;
        end
    end

    % search in the next level
    m = m - 1;
    v1 = treeSearch(v0 + vCurrent);
    m = m + 1;
    
    % update the bound
    minBound = max(minBound, v1);
    
    % unfill
    for ii = 1 : n
        ci = idxs{m, ii}(corr(m, ii));
        visG(ci, ii) = 1;
        
        for jj = ii + 1 : n
            cj = idxs{m, jj}(corr(m, jj));
            PGs(ci, cj, ii, jj) = 0;
        end
    end
    
    i = i + 1;
end
v = minBound;
