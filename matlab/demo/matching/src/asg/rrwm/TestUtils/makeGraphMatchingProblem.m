function problem = makeGraphMatchingProblem(settings, nControl)
%
% Input
%
% Output

% Control variable varification
if nargin < 2
    nOutlier = settings.nOutlier;
    ratioFill = settings.ratioFill;
    deformation = settings.deformation;
else
    eval([settings.Con.name ' = ' num2str(settings.Con.var(nControl)) ';']);
    for i = 1 : length(settings.Fix)
        eval([settings.Fix(i).name ' = ' num2str(settings.Fix(i).var) ';']);
    end
end

scaleSigma = settings.scaleSigma;

% Generate two graps
% define the number of nodes in both graph considering outliers
n1In = settings.nInlier;
n2In = settings.nInlier;
if settings.bOutlierOneSide
    n1Out = 0;
else
    n1Out = nOutlier;
end
n2Out = nOutlier;
n1 = n1In + n1Out;
n2 = n2In + n2Out; % n1 <= n2

% Generate graph 1 
G1 = tril(rand(n1), -1); % lower triangular graph
G1 = G1 + G1';
P = tril(rand(n1), -1);
P = P + P';
P = P > ratioFill;
G1(P) = NaN;

if settings.bPermuteGraph
    randIdx = randperm(n2); % 1 : n1In randIdx(1:n1In) will be true match
else
    randIdx = 1 : n2;
end

% deformation
N = deformation * tril(randn(n2), -1);
N = N + N';

% graph 2
G2 = tril(rand(n2), -1);
G2 = G2 + G2';
P = tril(rand(n2),-1);
P = P + P';
P = P > ratioFill;
G2(P) = NaN;
G2(randIdx(1 : n1In), randIdx(1 : n1In)) = G1(1 : n1In, 1 : n1In);
G2 = G2 + N;

E12 = ones(n1, n2);
[L12(:, 1) L12(:, 2)] = find(E12);

% affinity
M = zeros(n1 * n2);
for j = 1 : length(M) - 1
    for i = j + 1 : length(M)
        if isnan(G1(L12(j, 1), L12(i, 1))) || isnan(G2(L12(j, 2), L12(i, 2)))
            M(j, i) = 0;
        elseif ~settings.bUseConflict && xor(L12(j, 1) == L12(i, 1), L12(j, 2) == L12(i, 2))
            M(j, i) = 0;
        else
            M(j, i) = exp(-(G1(L12(j, 1), L12(i, 1)) - G2(L12(j, 2), L12(i, 2))) ^ 2 / scaleSigma);
        end
    end
end

target = accumarray([(1 : n1In)' randIdx(1 : n1In)'], 1, [n1 n2]);

% Return results
problem.G1 = G1;
problem.G2 = G2;
problem.E12 = E12;
problem.M = M + M';
problem.T12 = target;
