close all;
clear all;
clc
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% Get all the directory not yet post-processed
simuDirAll = '/calculSSD/salome/Simulation-10mai'; 

% simuDir = dir(sprintf('%s/simulation_rms_*', simuDirAll));
simuDir = dir(sprintf('%s/Bone*', simuDirAll));

dirFlags = [simuDir.isdir]; 
simuDir = simuDir(dirFlags); 

errorBis = [];
for idx = 1:numel(simuDir)
    simu_dir = simuDir(idx);
    if ~exist(fullfile(simu_dir.folder, simu_dir.name, 'postProcess.mat'))
        errorBis(end+1) = idx;
    end
end

%% Process possible directory, and save the index of those for which the simulation is not yet finished
errorRemain = [];
for idx = errorBis
    try
        simu_dir = simuDir(idx);
        [SpecularModel, SpecularProbaMap, OrientationtMap, reconstruction] = getPostProcess(fullfile(simu_dir.folder, simu_dir.name));
        save(fullfile(simu_dir.folder, simu_dir.name, 'postProcess.mat'), 'SpecularProbaMap', 'OrientationtMap', 'reconstruction');
        sprintf('Specular map computed for %s'  simu_dir.name);
    catch
        errorRemain(end+1) = idx;
    end
    
end

%% Post process with the same SpecularModel 
errorRemain = [];

% Get parameters of the simulation 
parameters = load(fullfile(simuDir(1).folder, simuDir(1).name, 'parameters.mat'));
recorded = LoadRfData(parameters.probe, fullfile(simuDir(1).folder, simuDir(1).name));

% Compute parameters required to reconstruct the image using DAS and/or specular transform
[acquisition, reconstruction] = GenerateParamRecon(recorded, parameters);

% Compute the specular transform
[SpecularTransform, TiltAngles] = ComputeSpecularTransform(reconstruction, acquisition, parameters);
%%
for idx = errorBis
    try
        simu_dir = simuDir(idx);
        plotMap = false;
        % Get parameters of the simulation 
        parameters = load(fullfile(simu_dir.folder, simu_dir.name, 'parameters.mat'));
        recorded = LoadRfData(parameters.probe, fullfile(simu_dir.folder, simu_dir.name));
        
        % Compute parameters required to reconstruct the image using DAS and/or specular transform
        [acquisition, reconstruction] = GenerateParamRecon(recorded, parameters);
        [SpecularModel, SpecularProbaMap, OrientationtMap] = ComputeSpecularityModel(SpecularTransform, ...
            acquisition, reconstruction, TiltAngles, plotMap, parameters, fullfile(simu_dir.folder, simu_dir.name));
        save(fullfile(simu_dir.folder, simu_dir.name, 'postProcess.mat'), 'SpecularProbaMap', 'OrientationtMap', 'reconstruction');
    catch
        errorRemain(end+1) = idx;
    end
    
end

%%
simuDirAll = '/calculSSD/salome/Simulation-10mai'; 
% simuDirAll = '~/Documents/BoneRugosity/SIMSONIC/Simulation/Simulation-19avr';
% simuDir = dir(sprintf('%s/simulation_rms_*', simuDirAll));
simuDir = dir(sprintf('%s/Bone227G*', simuDirAll));

dirFlags = [simuDir.isdir]; 
simuDir = simuDir(dirFlags); 


for idx = 1:numel(simuDir)
    simu_dir = simuDir(idx);
    if exist(fullfile(simu_dir.folder, simu_dir.name, 'postProcess.mat'))
        
        parameters = load(fullfile(simu_dir.folder, simu_dir.name, 'parameters.mat'));
        recorded = LoadRfData(parameters.probe, fullfile(simu_dir.folder, simu_dir.name));
        save(fullfile(simu_dir.folder, simu_dir.name, 'recorded.mat'), 'recorded')
    end
end



