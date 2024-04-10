clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% PROCESSING OF THE RF DATA 
% Generate a table to stock the geometry and the specularity Map
pathName = '/calculSSD/salome/Simulation-04avr';
rmsAll = 0.03 + (0:9) * 0.05;      % List of rms values (mm)
corrAll = [0.5 1 2 4];             % List of correlation length (mm)
MapAll = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));
SpecuAll = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));

%% PLOT ALL GEOMETRY MAP
for i = 1:numel(rmsAll)
    rms = str2double(MapAll.Properties.RowNames{i});
    for j = 1:numel(corrAll)
        % Generate path to simulation directory
        corr = str2double(MapAll.Properties.VariableNames{j});
        format = 'simulation_rms_%.2f_cl_%.1f/';
        simu_dir = [pathName, sprintf(format, rms, corr)];

        % Recover simulation Map
        [Map,~,~] = SimSonic2DReadMap2D([simu_dir, 'Geometry.map2D']);
        MapAll{i, j} = {Map};
    end
end

%% PLOT ALL MAPS
figure
t = tiledlayout(numel(rmsAll),numel(corrAll),'TileSpacing','tight');
% Plot Map
for i = 1:numel(rmsAll)
    for j = 1:numel(corrAll)
        % Generate path to simulation directory
        nexttile(t)
        img = MapAll{i,j}{1}(600:1200, :);
        imagesc(img);
%         axis equal
        % Remove x and y axes
        set(gca, 'xtick', [], 'ytick', []);
        
        if i == 1
            title(MapAll.Properties.VariableNames{j}, 'FontSize', 8, 'FontWeight', 'normal');
        end
        if j == 1
            ylabel(MapAll.Properties.RowNames{i}, 'FontSize', 8);
        end
    end
end
title(t, 'Correlation length (mm)', 'FontSize', 10)
ylabel(t, 'Root mean square (mm)', 'FontSize', 10)
xlabel(t, 'Interface profile for various RMS and correlation length', 'FontSize', 15, 'FontWeight', 'bold')

%% PLOT ALL SPECULARITY MAPS 
% Collect all simulations directory for the type of test.
SpecuAll = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));

for simul_dir = pathAllSimulation

end
