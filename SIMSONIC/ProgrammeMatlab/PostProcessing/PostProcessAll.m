close all;
clear all;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% COMPUTE SPECULARITY MAP FOR ALL SIMULATION
simuDirAll = '/calculSSD/salome/Simulation-04avr'; 
simuDir = dir(sprintf('%s/simulation_rms_*', simuDirAll));

dirFlags = [simuDir.isdir]; 
simuDir = simuDir(dirFlags); 

errorSimu = [];
for idx = 23:numel(simuDir)
    try
        simu_dir = simuDir(idx);
        [SpecularModel, SpecularProbaMap, OrientationtMap, reconstruction] = PostProcessing(fullfile(simu_dir.folder, simu_dir.name));
        save(fullfile(simu_dir.folder, simu_dir.name, 'postProcess.mat'), 'SpecularModel', 'SpecularProbaMap', 'OrientationtMap', 'reconstruction');
    catch
        errorSimu(end+1) = idx;
    end
    
end

