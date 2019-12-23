%% Implementaion of Leordeanu Marius' Spectral Matching method
% by Jungmin Lee in CVL
function X = SpectralMatching(M)

% Get eigenVector with largest eigenValue
options.disp = 0; 
M = sparse(M);
[eigenVectorMat eigenValue] = eigs(M, 1, 'lm', options); % EigenSolver
eigenVector = eigenVectorMat(:, 1); % Principal eigenvector
X = abs(eigenVector);
