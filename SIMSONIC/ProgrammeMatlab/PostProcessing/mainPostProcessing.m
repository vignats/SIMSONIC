clear; 
close all;
clc;
addpath(genpath('~/Documents')); 

%% PROCESSING OF THE RF DATA 
% Directory of the simulation, uncomment the needed information regarding
% the type of simulation
pathname = '~/Documents/BoneRugosity/SIMSONIC/Simulation/';

% Gaussian distribution of rough interface
rms = 0.5;
corr = 0.1;
format = 'simulation_rms_%.1f_cl_%.1f/';
simulation_name = sprintf(format, rms, corr);
simu_dir = [pathname, simulation_name];

% Ex-vivo rough interface
bone = '227G';
image = '1590';
fc = 1.25;
simulation_name = ['Bone', bone, '-Image', image, '-F', num2str(fc), '/'];
simu_dir = [pathname, simulation_name]; 

% Get parameters of the simulation 
parameters = load(fullfile(simu_dir, 'parameters.mat'));
recorded = LoadRfData(parameters.probe, simu_dir);

% Get Map and interface of the corresponding simulation
profile = load(fullfile(simu_dir, 'interface.mat'));

%% Compute parameters required to reconstruct the image using DAS and/or specular transform
[acquisition, reconstruction] = GenerateParamRecon(recorded);

%% Compute the specular transform
[SpecularTransform, TiltAngles] = ComputeSpecularTransform(reconstruction, acquisition);

%% Compute teh specularity map and angle
plotMap = true;
[SpecularModel, OrientationMap, SpecularityMap] = ComputeSpecularityModel(SpecularTransform, reconstruction, acquisition, TiltAngles, simu_dir, plotMap);

%% Plot map
       
figure('Position',[1317 1 1244 1321] ),
[Map,Nx,Nz] = SimSonic2DReadMap2D([simu_dir, 'Geometry.map2D']);
% nexttile(tlo,2)
X = 5:0.01:25-0.01; X = X -mean(X);
Z = flipud(4:0.01:15-0.01);
imagesc(X, Z, Map(400:end-1, 500:2499))
axis equal
xlabel('Width (mm)');
ylabel('Depth (mm)');
title(sprintf('Simulation map, rms = %f, corr = %f', parameters.interface.rms, parameters.interface.corr));
