clear; 
close all;
clc;
addpath(genpath('~/Documents')); 

%% PROCESSING OF THE RF DATA 
% Generate a table to stock the Map
rmsAll = 0.03 + (0:9) * 0.05;      % List of rms values (mm)
corrAll = [0.5 1 2 4];             % List of correlation length (mm)
MapAll = table('Size', [numel(corrAll), numel(rmsAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(rmsAll)), ...
                'VariableNames', cellstr(string(rmsAll)), ...
                'RowNames', cellstr(string(corrAll)));

figure
t = tiledlayout(numel(corrAll),numel(rmsAll),'TileSpacing','Compact');
for i = 1:numel(corrAll)
    corr = str2double(MapAll.Properties.RowNames{i});
    for j = 1:numel(rmsAll)
        % Generate path to simulation directory
        rms = str2double(MapAll.Properties.VariableNames{j});
        pathName = '~/Documents/BoneRugosity/SIMSONIC/Simulation/';
        format = 'simulation_rms_%.2f_cl_%.1f/';
        simu_dir = [pathName, sprintf(format, rms, corr)];

        % Recover simulation Map
        [Map,~,~] = SimSonic2DReadMap2D([simu_dir, 'Geometry.map2D']);
        MapAll{i, j} = {Map};
        % Plot Map
        nexttile
        plot(MapAll{i,j}{1}); % Vous devez personnaliser cette partie en fonction de vos données

        title(sprintf('corr = %s & rms = %s', MapAll.Properties.RowNames{i}, MapAll.Properties.VariableNames{j}));    
    end
end

%%
% Collect all simulations directory for the type of test.
pathSimulation = '~/Documents/BoneRugosity/SIMSONIC/Simulation/';
pathAllSimulation = dir(sprintf('%s/simulation_rms_*', pathSimulation));
% Ensure that DirAllSimulation only contain directories. 
dirFlags = [pathAllSimulation.isdir]; 
pathAllSimulation = pathAllSimulation(dirFlags); 

for simul_dir = pathAllSimulation
    
	% Get parameters of the simulation 
    parameters = load(fullfile(simu_dir, 'parameters.mat'));
    recorded = LoadRfData(parameters.probe, simu_dir);
    
    % Get Map and interface of the corresponding simulation
    profile = load(fullfile(simu_dir, 'interface.mat'));
    
    % Compute parameters required to reconstruct the image using DAS and/or specular transform
    [acquisition, reconstruction] = GenerateParamRecon(recorded);
    
    % Compute the specular transform
    [SpecularTransform, TiltAngles] = ComputeSpecularTransform(reconstruction, acquisition);
    
    % Compute teh specularity map and angle
    plotMap = false;
    [SpecularModel, OrientationMap, SpecularityMap] = ComputeSpecularityModel(SpecularTransform, reconstruction, acquisition, TiltAngles, simu_dir, plotMap);
    
    % Save the needed Map in the file
    if saveParam
        save(fullfile(simu_dir, 'SpecularityMap.mat'), OrientationMap, SpecularityMap)
    end
end
