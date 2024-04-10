close all;
clear all;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% COMPUTE SPECULARITY MAP FOR ALL SIMULATION
simuDirAll = '/calculSSD/salome/Simulation-04avr'; 
simuDir = dir(sprintf('%s/simulation_rms_*', simuDirAll));

dirFlags = [simuDir.isdir]; 
simuDir = simuDir(dirFlags); 

% errorSimu = [];
errorBis = [];
for idx = errorSimu
    try
        simu_dir = simuDir(9);
        [SpecularModel, SpecularProbaMap, OrientationtMap, reconstruction] = PostProcessing(fullfile(simu_dir.folder, simu_dir.name));
        save(fullfile(simu_dir.folder, simu_dir.name, 'postProcess.mat'), 'SpecularProbaMap', 'OrientationtMap', 'reconstruction');
    catch
        errorBis(end+1) = idx;
    end
    
end

