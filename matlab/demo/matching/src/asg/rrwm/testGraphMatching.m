clear variable;
close all;
clc;
initPath;

%% Problem settings
settings.nInlier = 15;

switch 1
    case 1 % varing deformation
        settings.nOutlier = 0;
        settings.deformation = 0 : 0.05 : 0.4; % varying deformation
        settings.ratioFill = 1; 
    case 2 % varing outlier
        settings.nOutlier = 0 : 2 : 30; % varying # of outliers
        settings.deformation = 0;
        settings.ratioFill = 1;
    case 3 % varing edge density
        settings.nOutlier = 10;
        settings.deformation = 0.1;
        settings.ratioFill = 0.3 : 0.1 : 1; % varying edge density
    otherwise
end

settings.nTest = 10; % # of test trials
settings.scaleSigma = 0.15;
settings.bOutlierOneSide = 1; % 0: both 1: domain B
settings.bPermuteGraph = 1;
settings.bNormalizeAll = 0; % 0: method dependent  1: for all methods  -1: for no method
settings.bUseConflict = 0; % 0: delete conficting element  1: use conflicting elements

settings.bMulti = 0; % 0: single  1: double
settings.bVector = 0; % 0: distance  1: vector

settings.bDisplayLegend = 0;

%% set configuration for algorithms
configMethod;
setVariable;

%% Test part
score = zeros(length(settings.Con.var), length(settings.Algorithm), settings.nTest);
time = zeros(length(settings.Con.var), length(settings.Algorithm), settings.nTest);
energy = zeros(length(settings.Con.var), length(settings.Algorithm), settings.nTest);
energyIn = zeros(length(settings.Con.var), length(settings.Algorithm), settings.nTest);

for i = 1 : settings.nTest
    fprintf('Test %d of %d ', i, settings.nTest);
    for j = 1 : length(settings.Con.var)
        problem = makeGraphMatchingProblem(settings, j); % please input control number!!
        [tempScore tempEnergy tempInEnergy tempTime] = solveProblem(problem, settings);
        for k = 1:length(settings.Algorithm)
            score(j, k, i) = tempScore(k);
            energy(j, k, i) = tempEnergy(k);
            energyIn(j, k, i) = tempInEnergy(k);
            time(j, k, i) = tempTime(k);
        end
        fprintf('.');
    end
    fprintf('\n');
end

meanScore = mean(score, 3) / settings.nInlier;
meanEnergy = mean(energy, 3);
meanInEnergy = mean(energyIn, 3);
meanTime = mean(time, 3);

%% Show and Save
show;
