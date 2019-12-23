function [nTrueMatch objValue objInValue computeTime] = solveProblem(problem, settings)
% Solving the matching problem
% Jungmin Lee in CVL
% problem : structure that contains all information

Algorithm = settings.Algorithm;

% Solve the problem using various methods
nTrueMatch = zeros(length(Algorithm), 1);
objValue = zeros(length(Algorithm), 1);
objInValue = zeros(length(Algorithm), 1);
computeTime = zeros(length(Algorithm), 1);

for k = 1 : length(Algorithm)
    tic;
    X = feval(Algorithm(k).fhandle, problem, Algorithm(k));
    computeTime(k) = toc;
    [nTrueMatch(k) objValue(k) objInValue(k)] = checkAnswer(X, problem.M, problem.T12, Algorithm(k));
end
