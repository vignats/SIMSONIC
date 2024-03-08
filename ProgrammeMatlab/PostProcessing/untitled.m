clear; 
close all;
clc;
addpath(genpath('~/Documents')); 

%% PROCESSING OF THE RF DATA 
% Directory of the simulation
rms = 0.5;
corr = 0.1;

filename = '~/Documents/SIMSONIC/Simulation/';
format = 'simulation_rms_%.1f_cl_%.1f/';
simulation_name = 'simulation_rough_interface/'; % sprintf(format, rms, corr);
simu_dir = [filename, simulation_name];

% Get parameters of the simulation 
parameters = load(fullfile(simu_dir, 'parameters.mat'));
recorded = LoadRfData(parameters.probe, simu_dir);

% Get Map of the corresponding simulation
[Map,Nx,Ny] = SimSonic2DReadMap2D([simu_dir, 'Geometry.map2D']);

%% COMPUTE SPECULAR TRANSFORM
% 1. Compute the delayed signal 
% 2. Plot the DAS image
% 3. Compute the specular transform for all pixel
% 4. Plot it for a pixel of interest (x_pixel, z_pixel)
% 5. Compute the model of specular transform
% 6. Plot the specularity map and the original map

pixel.z = 8e-3;
pixel.x = 0e-3;
[FULL_DELAYED_MAP_OF_RF_DATA_PIXEL, DAS_IMAGE_dB, SPECULAR_TRANSFORM] = ReconstructImage(simu_dir, recorded, Map, parameters.grid, pixel);
