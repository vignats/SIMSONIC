% CREATE NECESSARY FILE FOR A SIMULATION 
clear; 
close all;
clc;
addpath(genpath('~/Documents')); 

%% ACQUISITION PARAMETERS 
% If the saveParam argument is true, a
% 'parameters.mat' file will be registred in the directory.
saveParam = false;
[param, grid, probe, medium, interface, signal, simu_dir] = GenerateAllParameters(saveParam);

%% COMPUTATION OF THE SIGNAL AND GEOMETRY FILES
fprintf('---------- Geometry and Signal genration ----------\n');
print_plot = true;               % To plot the signal and the geometry.
print_hist = false;              % To plot the size and position distributions of the pore

% If the path to the simulation file is saved as an argument, a
% the signal and geometry files will be registred in the directory.
MakeSgl(param, grid, medium, signal, print_plot, simu_dir)
% MakeGeometryRugosity(grid, probe, medium, interface, print_plot, simu_dir);
MakeGeometryInterface(grid, probe, medium, interface, print_plot);

%% CREATION OF THE PARAMETERS FILE IN EACH TRANSMITTER DIRECTORY
fprintf('---------- Parameters genration ----------\n');
print_plot = false;
MakeParameters(param, grid, probe, medium, simu_dir, print_plot)