function problem = makePointMatchingProblem(settings, nControl)
% Control variable varification

if nargin < 2
    nOutlier = settings.nOutlier;
    deformation = settings.deformation;
else
    eval([settings.Con.name ' = ' num2str(settings.Con.var(nControl)) ';']);
    for i = 1:length(settings.Fix)
        eval([settings.Fix(i).name ' = ' num2str(settings.Fix(i).var) ';']);
    end
end

% Generate two graphs
nInlier = settings.nInlier;
sizeRegion = 256 * sqrt((nInlier+nOutlier)/10);

if ~settings.bOutlierOneSide
    n1 = nInlier + nOutlier; % number of node in Graph 1
    n2 = nInlier + nOutlier; % number of node in Graph 2
    node1 = sizeRegion * rand(n1, 2); % Random Points in 2D space
    node2 = node1(1 : nInlier, :) + random('normal', 0, deformation, nInlier, 2);
    node2 = [node2; sizeRegion * rand(nOutlier, 2)];
    
elseif settings.bOutlierOneSide
    n1 = nInlier; % number of node in Graph 1
    n2 = nInlier + nOutlier; % number of node in Graph 2
    node1 = sizeRegion * rand(n1, 2); % Random Points in 2D space
    node2 = node1(1 : nInlier, :) + random('normal', 0, deformation, nInlier, 2);
    node2 = [node2; sizeRegion * rand(nOutlier, 2)];
end

%% Permute Graph2
if settings.bPermuteGraph
    permIdx = randperm(n2);
else
    permIdx = 1:n2;
end
node2(permIdx,:) = node2;

%% Construct matching indicator matrix E12
E12 = ones(n1,n2);
n12 = nnz(E12);
% [L12(:,2) L12(:,1)] = find(E12);
[L12(:,1) L12(:,2)] = find(E12);

%% Construct Affinity matrix M - using distance
E1 = ones(n1); E2 = ones(n2);
[L1(:,1) L1(:,2)] = find(E1);
[L2(:,1) L2(:,2)] = find(E2);
G1 = node1(L1(:, 1), :) - node1(L1(:,2),:);
G2 = node2(L2(:, 1), :) - node2(L2(:,2),:);
G1 = sqrt(G1(:,1).^2+G1(:,2).^2);
G2 = sqrt(G2(:,1).^2+G2(:,2).^2);
G1 = reshape(G1, n1, n1);
G2 = reshape(G2, n2, n2);

M = zeros(n12);
for j = 1:n12-1
    for i = j+1:n12
        dist1 = norm(node1(L12(j,1),:)-node1(L12(i,1),:));
        dist2 = norm(node2(L12(j,2),:)-node2(L12(i,2),:));
        M(j,i) = exp(-(dist1-dist2)^2/settings.scaleSigma); 
        M(i,j) = M(j,i);
    end
end

% True matches
T12 = zeros(n1,n2);
for i = 1:nInlier
    T12(i,permIdx(i)) = 1;
end
node = [node1; node2];
sizeRegion = [min(node(:,1)) max(node(:,1)) min(node(:,2)) max(node(:,2))];

% Return results
problem.sizeRegion = sizeRegion;
problem.node1 = node1;
problem.node2 = node2;
problem.G1 = G1;
problem.G2 = G2;
problem.scaleSigma = settings.scaleSigma;
problem.E12 = E12;
problem.M = M;
problem.T12 = T12;
