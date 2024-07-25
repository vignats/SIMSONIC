% CREATE NECESSARY FILE FOR A SIMULATION 
clear; 
close all;
clc;
addpath(genpath('~/Documents')); 
addpath(genpath('/calculSSD/salome'));

%% ACQUISITION PARAMETERS 
% If the saveParam argument is true, a
% 'parameters.mat' file will be registred in the directory.
saveParam = false;
[param, grid, probe, medium, interface, signal, simu_dir] = GenerateAllParameters(saveParam);

%% COMPUTATION OF THE SIGNAL AND GEOMETRY FILES
print_hist = true;              % To plot the size and position distributions of the pore
verify = 'plot';
% If the path to the simulation file is saved as an argument, a
% the signal and geometry files will be registred in the directory.
% MakeSgl(param, grid, medium, signal, false, simu_dir)
Map = MakeGeometryRugosity(grid, probe, medium, interface, true, print_hist);
% MakeGeometry(grid, probe, medium, interface, verify);

%% CREATION OF THE PARAMETERS FILE IN EACH TRANSMITTER DIRECTORY
fprintf('---------- Parameters genration ----------\n');
print_plot = false;
MakeParameters(param, grid, probe, medium, simu_dir, print_plot)