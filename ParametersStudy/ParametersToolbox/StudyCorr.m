% This script aim to study the relevance of the ComputeInterfaceParameters function
% regarding the MakeGeometryInterface function
clear; 
close all;
clc;
addpath(genpath('~/Documents')); 
%% Interface Generation 
saveParam = false;
[param, grid, probe, medium, interface, signal, simu_dir] = GenerateAllParameters(saveParam);

print_plot = false;             % To plot the signal and the geometry

% Compute the correlation length with many repetition to study the confidence interval
nbRepetition = 1000;
parametersValues = zeros(3, nbRepetition);
corrValues = zeros(1, nbRepetition);

for i = 1:nbRepetition
    [Map, heights] = MakeGeometryInterface(grid, probe, medium, interface, print_plot);
    corrValues(1, i) = ComputeCorr(heights, grid.step);

    [Rq, Corr, rugosity] = ComputeInterfaceParameters(Map, grid, probe, medium);
    parametersValues(1, i) = Rq;
    parametersValues(2, i) = Corr;
    parametersValues(3, i) = rugosity;
    disp(i)
end 

figure;
histogram(corrValues);
xlabel('Correlation length (mm)');
ylabel('Number of occurences');
title('Computed correlation length distribution with profile', sprintf('Initial correlation length = %.02f mm', interface.corr));

figure;
histogram(parametersValues(2, :));
xlabel('Correlation length (mm)');
ylabel('Number of occurences');
title('Computed correlation length distribution with Map', sprintf('Initial correlation length = %.02f mm', interface.corr));

figure;
histogram(parametersValues(1, :));
xlabel('Root mean square (mm)');
ylabel('Number of occurences');
title('Root mean square with Map', sprintf('Initial rms = %.02f mm', interface.rms));

figure;
histogram(parametersValues(2, :));
xlabel('Rugosity (%)');
ylabel('Number of occurences');
title('Rugosity (%) distribution with Map');

stdParam = zeros(3,1); 
stdParam(1, 1) = std(parametersValues(1, :));
stdParam(2, 1) = std(parametersValues(2, :));
stdParam(3, 1) = std(parametersValues(3, :));
stdProf = std(corrValues);

meanParam = zeros(3,1); 
meanParam(1, 1) = mean(parametersValues(1, :));
meanParam(2, 1) = mean(parametersValues(2, :));
meanParam(3, 1) = mean(parametersValues(3, :));
meanProf = mean(corrValues);