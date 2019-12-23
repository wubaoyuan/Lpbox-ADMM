%% MATLAB demo script for point mathcing
% This is for demonstrating Reweighted Random Walks Graph Matching ECCV2010
% http://cv.snu.ac.kr/research/~RRWM/
%
% Minsu Cho, Jungmin Lee, and Kyoung Mu Lee, 
% Reweighted Random Walks for Graph Matching, 
% Proc. European Conference on Computer Vision (ECCV), 2010
%
% written by Minsu Cho & Jungmin Lee 2010, Seoul National University, Korea
% http://cv.snu.ac.kr/~minsucho/
% http://cv.snu.ac.kr/~jungminlee/

clear all; close all; clc; initPath;

%% Problem settings
settings.nInlier = 15;

switch 1
    case 1 % varying deformation
        settings.deformation = 0:1:20; % varying deformation 
        settings.nOutlier = 0;
    case 2 % varying num of outliers
        settings.deformation = 0; 
        settings.nOutlier = 0:2:30; % varying # of outliers
    otherwise
end

settings.nTest = 5; % # of test trials

settings.bOutlierOneSide = 0; % 0: both  1: only at domain B
settings.bPermuteGraph = 1;
settings.scaleSigma = 100;
settings.ratioFill = [];

%% set configuration for algorithms
configMethod;
setVariable;

%% Testing part
score = zeros(length(settings.Con.var), length(settings.Algorithm), settings.nTest);
time = zeros(length(settings.Con.var), length(settings.Algorithm), settings.nTest);
energy = zeros(length(settings.Con.var), length(settings.Algorithm), settings.nTest);

for i = 1:settings.nTest
    fprintf('Test %d of %d ', i, settings.nTest);
    for j = 1:length(settings.Con.var)
        problem = makePointMatchingProblem( settings, j); % please input control number!!
        [ tempScore tempEnergy tempGarbage tempTime ] = solveProblem( problem, settings );
        for k = 1:length(settings.Algorithm)
            score(j,k,i) = tempScore(k);
            energy(j,k,i) = tempEnergy(k);
            time(j,k,i) = tempTime(k);
        end
        fprintf('.');
    end
    fprintf('\n');
end

meanScore = mean(score, 3)/settings.nInlier;
meanEnergy = mean(energy, 3);
meanTime = mean(time, 3);

stdScore = std(score, 0, 3);
stdEnergy = std(energy, 0, 3);
stdTime = std(time, 0, 3);

clear problem tempScore tempEnergy tempGarbage tempTime

%% Show results
show; clear i j k