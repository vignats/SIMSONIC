clear; 
close all;
clc;
addpath(genpath('~/Documents')); 

%% PROCESSING OF THE RF DATA 
% Directory of the simulation, uncomment the needed information regarding
% the type of simulation
filename = '~/Documents/BoneRugosity/SIMSONIC/Simulation/';

% Gaussian distribution of rough interface
% rms = 0.5;
% corr = 0.1;
% format = 'simulation_rms_%.1f_cl_%.1f/';
% simulation_name = sprintf(format, rms, corr);
% simu_dir = [filename, simulation_name];

% Ex-vivo rough interface
bone = '227G';
image = '1590';
fc = 1.25;
simulation_name = ['Bone', bone, '-Image', image, '-F', num2str(fc), '/'];
simu_dir = [filename, simulation_name]; 

% Get parameters of the simulation 
parameters = load(fullfile(simu_dir, 'parameters.mat'));
recorded = LoadRfData(parameters.probe, simu_dir);

% Get Map and interface of the corresponding simulation
[Map,Nx,Ny] = SimSonic2DReadMap2D([simu_dir, 'Geometry.map2D']);
profile = load(fullfile(simu_dir, 'interface.mat'));

%% COMPUTE SPECULAR TRANSFORM
% 1. Compute the delayed signal 
% 2. Plot the DAS image
% 3. Compute the specular transform for all pixel
% 4. Plot it for a pixel of interest (x_pixel, z_pixel)
% 5. Compute the model of specular transform
% 6. Plot the specularity map and the original map

pixel.z = 8e-3;
pixel.x = 0e-3;
[FULL_DELAYED_MAP_OF_RF_DATA_PIXEL, DAS_IMAGE_dB, SPECULAR_TRANSFORM] = ReconstructImage(simu_dir, recorded, Map, profile, pixel);