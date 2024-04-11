close all;
clear all;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% COMPUTE SPECULARITY MAP FOR ALL SIMULATION
simuDirAll = '/calculSSD/salome/Simulation-04avr'; 
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

%%
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

