close all;
clear all;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% Get all the directory not yet post-processed
simuDirAll = '/calculSSD/salome/Simulation-19avr'; 
simuDir = dir(sprintf('%s/simulation_rms_*', simuDirAll));

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
        [SpecularModel, SpecularProbaMap, OrientationtMap, reconstruction] = PostProcessing(fullfile(simu_dir.folder, simu_dir.name));
        save(fullfile(simu_dir.folder, simu_dir.name, 'postProcess.mat'), 'SpecularProbaMap', 'OrientationtMap', 'reconstruction');
    catch
        errorRemain(end+1) = idx;
    end
    
end



