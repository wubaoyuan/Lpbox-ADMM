function X = SM(problem, AlgConfig)
%% The principal eigenvector of the affinity matrix
X = SpectralMatching(problem.M); 

%% Get index from principle eigenVector
X = feval(AlgConfig.postProcess, problem, X);
