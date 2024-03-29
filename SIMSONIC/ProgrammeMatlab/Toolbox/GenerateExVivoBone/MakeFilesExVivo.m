% CREATE NECESSARY FILE FOR A SIMULATION 
clear; 
close all;
clc;
addpath(genpath('~/Documents')); 

%% ACQUISITION PARAMETERS 
% If the path to the simulation file is saved as an argument, a
% 'parameters.mat' file will be registred in the directory.
saveParam = false;
[param, grid, probe, medium, interface, signal, filter, simu_dir] = GenerateAllParametersExVivo(false);

%% COMPUTATION OF THE SIGNAL AND GEOMETRY FILES
fprintf('---------- Geometry and Signal genration ----------\n');
print = true;            % To plot the signal and the geometry.
print_hist = true;       % To plot the size and position distributions of the pore

MakeSgl(param, grid, medium, signal, print, simu_dir);
MakeGeometryExVivo(grid, interface, filter, print);

%% CREATION OF THE PARAMETERS FILE IN EACH TRANSMITTER DIRECTORY
fprintf('---------- Parameters genration ----------\n');
print = false;
MakeParameters(param, grid, probe, medium, simu_dir, print)