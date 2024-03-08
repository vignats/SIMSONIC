% CREATE NECESSARY FILE FOR A SIMULATION 
clear; 
close all;
clc;
addpath(genpath('~/Documents')); 

%% CREATION OF THE SIMULATION DIRECTORY
filename = '~/Documents/SIMSONIC/Simulation/';
format = 'simulation_rms_%.1f_cl_%.1f/';
simulation_name = 'test/'; %sprintf(format, interface.rms, interface.cl);
simu_dir = [filename, simulation_name]; 
if ~exist(simu_dir,'dir')
		mkdir(simu_dir);
end

%% ACQUISITION PARAMETERS 
% If the path to the simulation file is saved as an argument, a
% 'parameters.mat' file will be registred in the directory.
GenerateAllParameters(simu_dir);

%% COMPUTATION OF THE SIGNAL AND GEOMETRY FILES
fprintf('---------- Geometry and Signal genration ----------\n');
print = true;            % To plot the signal and the geometry.
print_hist = true;       % To plot the size and position distributions of the pore

MakeSgl(param, grid, medium, signal, simu_dir, print)
% MakeGeometry(grid, probe, medium, interface, simu_dir, print, print_hist)
% Map_50 = MakeGeometryRugosity(grid, probe, medium, interface, simu_dir, print, print_hist);
Map_I = MakeGeometryInterface(grid, probe, medium, interface, simu_dir, print);

%% CREATION OF THE PARAMETERS FILE IN EACH TRANSMITTER DIRECTORY
fprintf('---------- Parameters genration ----------\n');
print = false;
MakeParameters(param, grid, probe, medium, simu_dir, print)