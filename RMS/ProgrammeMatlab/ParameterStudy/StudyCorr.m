% This script aim to study the relevance of the ComputeCorr function
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
corrValues = zeros(1, nbRepetition);

for i = 1:nbRepetition
    [~, heights] = MakeGeometryInterface(grid, probe, medium, interface, print_plot);
    corrValues(1, i) = ComputeCorr(heights, grid.step);
end 

figure;
histogram(corrValues);
xlabel('Correlation length (mm)');
ylabel('Number of occurences');
title('Computed correlation length distribution', sprintf('Initial correlation length = %.02f mm', interface.corr));