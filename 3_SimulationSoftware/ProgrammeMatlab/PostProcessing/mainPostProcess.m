clear; 
close all;
clc;
addpath(genpath('~/Documents')); 
addpath(genpath('/calculSSD/salome'));

%% PROCESSING OF THE RF DATA 
% Directory of the simulation, uncomment the needed information regarding
% the type of simulation
% Gaussian distribution of rough interface
pathname = '/calculSSD/salome/Simulation-17juin/';
Rms = 0.43;
corr = 2;
format = 'simulation_rms_%.2f_cl_%.1f/';
simulation_name = sprintf(format, Rms, corr);
simu_dir = [pathname, simulation_name];
%%
% fc = 1.25;
% simulation_name = ['Bone', bone, '-Image', image, '-F', num2str(fc), '/'];

% Ex-vivo rough interface
% pathname = '/calculSSD/salome/Simulation-10mai/';
% bones = {'245D', '227G', '267G'};
% slices = {{0000, 1134, 3195, 3852, 5511}, {0000, 2002, 3595, 5614, 6721}, ...
%     {0000, 1700, 3410, 5478, 6716}}; 
% 
% i = 1; j = 1;
% bone.id = bones{i};
% bone.image = slices{i}{j};
% 
% simulation_name = ['Bone', bone.id, '-Image', sprintf('%04d', bone.image), '/'];
% simu_dir = fullfile(pathname, simulation_name); 

% Get parameters of the simulation 
parameters = load(fullfile(simu_dir, 'parameters.mat'));
recorded = LoadRfData(parameters.probe, simu_dir);

%% Compute parameters required to reconstruct the image using DAS and/or specular transform
[acquisition, reconstruction] = GenerateParamRecon(recorded, parameters);

%% Compute the specular transform
[SpecularTransform, TiltAngles] = ComputeSpecularTransform(reconstruction, acquisition, parameters);

%% Compute the specularity map and angle
plotMap = true;
[SpecularModel] = ComputeSpecularityModel(SpecularTransform, acquisition, reconstruction, TiltAngles, plotMap, parameters, simu_dir);

%% Plot map
       
figure('Position',[1317 1 1244 1321] ),
[Map,Nx,Nz] = SimSonic2DReadMap2D(fullfile(simu_dir, 'Geometry.map2D'));
% nexttile(tlo,2)
X = 5:0.01:25-0.01; X = X -mean(X);
Z = flipud(4:0.01:15-0.01);
imagesc(X, Z-4, Map(200:end-1, 500:2499))
axis image
xlabel('Width (mm)');
ylabel('Depth (mm)');
title('Simulation map')%, sprintf('rms = %.2f, corr = %.1f', parameters.interface.rms, parameters.interface.corr));
