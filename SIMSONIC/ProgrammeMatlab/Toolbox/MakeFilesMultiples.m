% CREATE NECESSARY FILE FOR A SIMULATION 
clear; 
close all;
clc;
addpath(genpath('~/Documents')); 
addpath(genpath('/calculSSD/salome'));

%% PARAMETERS TO MODIFY AT EACH SIMULATION 
% ACQUISITION PARAMETERS 
% If the saveParam argument is true, a
% 'parameters.mat' file will be registred in the directory.

saveParam = true;
[param, grid, probe, medium, interface, signal, simu_dir, corrAll, rmsAll] = GenerateAllParametersMultiples(saveParam);

%% COMPUTATION OF THE SIGNAL AND GEOMETRY FILES
fprintf('---------- Geometry and Signal genration ----------\n');
verify = 'null';               % To plot the signal and the geometry, and or confirm the heigth profile.
                              % Either 'plot', 'confirm' or 'null'
print_hist = false;           % To plot the size and position distributions of the pore

for i = 1:length(corrAll)
    for j = 1:length(rmsAll)
        interface.corr = corrAll(i);          % Correlation length (mm)
        interface.rms = rmsAll(j);            % Rms height (mm)
        
        simuDirAll = '~/Documents/BoneRugosity/SIMSONIC/Simulation/Simulation-19avr';
        format = 'simulation_rms_%.2f_cl_%.1f/';
        simulation_name = sprintf(format, interface.rms, interface.corr);
        simu_dir = fullfile(simuDirAll, simulation_name); 

        % If the path to the simulation file is saved as an argument, a
        % the signal and geometry files will be registred in the directory.
        MakeSgl(param, grid, medium, signal, false, simu_dir)

        MakeGeometryInterface(grid, probe, medium, interface, verify, simu_dir);

        MakeParameters(param, grid, probe, medium, simu_dir, false)
    end
end