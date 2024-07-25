% CREATE NECESSARY FILE FOR A SIMULATION 
clear; 
close all;
clc;
addpath(genpath('~/Documents')); 
addpath(genpath('/calculSSD/salome')); 

%% ACQUISITION PARAMETERS 
% If the saveParam is true, a 'parameters.mat' file will be registred 
% in the directory.
[param, grid, probe, medium, signal] = GenerateAllParametersExVivo();

%% Computation of the files for all bones images   
bones = {'245D', '227G', '267G'};
slices = {{1134, 3195, 3852, 5511}, {2002, 3595, 5614, 6721}, {1700, 3410, 5478, 6716}}; 

% for i = 1:numel(bones)
%     bone.id = bones{i};         % Bone from ex-vivo files
i = 3;  
bone.id = bones{i};         % Bone from ex-vivo files
for j = 1:numel(slices{i})
    bone.image = slices{i}{j};        % Slice selected
    
    % CREATE FILE
    withoutPores = true;
    dirname = '/calculSSD/salome/Simulation-10mai';

    simulation_name = ['Bone', bone.id, '-Image', sprintf('%04d', bone.image), '/'];
    
    if withoutPores
        simulation_name = ['Bone', bone.id, '-Image0000/'];
    end
    simu_dir = fullfile(dirname, simulation_name); 
    
    if ~exist(simu_dir,'dir')
        mkdir(simu_dir);
    end

    % COMPUTATION OF THE SIGNAL AND GEOMETRY FILES 
    fprintf('---------- Geometry and Signal genration ----------\n');
    printPlot = true;            % To plot the signal and the geometry.
    rotate = false;               % If the bone needs to be rotated
    % If the path to the simulation file is saved as an argument, a
    % the signal and geometry files will be registred in the directory.
    MakeSgl(param, grid, medium, signal, false, simu_dir);
    [Map, grid] = MakeGeometryExVivo(grid, probe, bone, rotate, withoutPores, printPlot, simu_dir);

    % CREATION OF THE PARAMETERS FILE IN EACH TRANSMITTER DIRECTORY
    fprintf('---------- Parameters genration ----------\n');
    printParam = false;
    param.length = sqrt((grid.depth - probe.depth)^2 + grid.width^2)*2e3/max(medium.cp);
    
    MakeParameters(param, grid, probe, medium, simu_dir, printParam)
    
    % SAVE PARAMETERS AND USEFUL FUNCTION  
    saveParam = true;
    if saveParam
        save(fullfile(simu_dir, 'parameters.mat'), 'param', 'grid',...
            'probe', 'medium', 'signal', 'bone', 'simu_dir')
    end
end