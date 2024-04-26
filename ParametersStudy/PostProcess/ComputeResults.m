clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% PROCESSING OF THE RF DATA 
% Generate a table to stock the geometry and the specularity Map
pathSimuAll = '/calculSSD/salome/Simulation-19avr';
rmsAll = 0.03 + (0:9) * 0.05;      % List of rms values (mm)
corrAll = [0.5 1 2 4];             % List of correlation length (mm)

Results = struct();

%% Compute maps and statistics
Results.Map = ComputeMap('Map', rmsAll, corrAll, pathSimuAll, Results);
Results.SpecuMap = ComputeMap('SpecuMap', rmsAll, corrAll, pathSimuAll, Results);

%% Plot maps
MapTitle = 'Interface profile for various RMS and correlation length';
PlotAll('Map', Results, rmsAll, corrAll, pathSimuAll, MapTitle);

MapTitle = 'Specular probability for various RMS and correlation length';
PlotAll('SpecuMap', Results, rmsAll, corrAll, pathSimuAll, MapTitle);

%% Plot statistics
SpecuProba = ComputeMap('Stat', rmsAll, corrAll, pathSimuAll, Results);
Results.Stat = struct();
[Results.Stat.meanROI, Results.Stat.stdROI] = deal(table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'double'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll))));
Results.Stat.linearROI = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));

for i = 1:numel(rmsAll)
    for j = 1:numel(corrAll)
        Results.Stat.meanROI{i,j} = SpecuProba{i,j}{1}.meanROI;
        Results.Stat.stdROI{i,j} = SpecuProba{i,j}{1}.stdROI;
        % Results.Stat.linearROI{i,j} = {SpecuProba{i,j}{1}.linearROI};
    end
end

figure;
surf(corrAll, rmsAll, table2array(Results.Stat.meanROI), 'EdgeColor', 'none');
xlabel('Correlation length (mm)');
ylabel('Root mean square (mm)');
zlabel('Mean Specular Probability');
title('Mean Specular Probability in the Region of Interest');
colorbar

%%
MapTitle = 'Specular probability along the lateral position';
PlotAll('Stat', Results, rmsAll, corrAll, pathSimuAll, MapTitle);

%% Compute initial parameters 
Param = ComputeMap('InterfaceParam', rmsAll, corrAll, pathSimuAll, Results);
Results.InitParam = struct();
Results.InitParam.Corr = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll))); 
Results.InitParam.Rms = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll))); 

for i = 1:numel(rmsAll)
    for j = 1:numel(corrAll)
        Results.InitParam.Rms{i,j} = Param{i,j}{1}(1);
        Results.InitParam.Corr{i,j} = Param{i,j}{1}(2);
    end
end

%% Save results 
save(fullfile(pathSimuAll, 'Results.mat'), 'Results')
%% Useful functions
function[Table] = ComputeMap(MapType, rmsAll, corrAll, pathSimuAll, Results)
    % Compute parameters once 
    format = 'simulation_rms_%.2f_cl_%.1f/';
    simuDir1 = fullfile(pathSimuAll, sprintf(format, rmsAll(1), corrAll(1)));
    parameters = load(fullfile(simuDir1, 'parameters.mat'));

    Table = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));

    % Generate Map and probability regarding the case 
    for i = 1:numel(rmsAll)
        rms = str2double(Table.Properties.RowNames{i});
        for j = 1:numel(corrAll)
            corr = str2double(Table.Properties.VariableNames{j});
            simu_dir = fullfile(pathSimuAll, sprintf(format, rms, corr));
    
            % Compute the required type of simulation 
            switch MapType
                case 'Map'
                    [Map,~,~] = SimSonic2DReadMap2D([simu_dir, 'Geometry.map2D']);
                    Table{i,j}{1} = Map;
                   
                case 'SpecuMap'
                    if exist(fullfile(simu_dir, 'postProcess.mat'))
                        postProcess = load(fullfile(simu_dir, 'postProcess.mat'));
                        Table{i,j}{1} = postProcess.SpecularProbaMap;
                    end
                case 'Stat'
                    recorded = LoadRfData(parameters.probe, simuDir1);
                    [~, reconstruction] = GenerateParamRecon(recorded);
                    probability = ProbaROI(Results.SpecuMap{i,j}{1}, reconstruction, parameters, false);
                    Table{i,j}{1} = probability;
                case 'InterfaceParam'
                    [Rq, Corr, rugosity] = ComputeInterfaceParameters(Results.Map{i,j}{1}, parameters.grid, parameters.probe, parameters.medium);
                    Table{i,j}{1} = {Rq, Corr, rugosity};
                otherwise 
                    warning('Unexpected Map type. No data computed.')
            end
        end
    end
end

function[] = PlotAll(MapType, Struct, rmsAll, corrAll, pathSimuAll, MapTitle)
    Table = getfield(Struct, MapType);

    % Compute parameters once 
    format = 'simulation_rms_%.2f_cl_%.1f/';
    simuDir1 = fullfile(pathSimuAll, sprintf(format, rmsAll(1), corrAll(1)));
    parameters = load(fullfile(simuDir1, 'parameters.mat'));

    figure
    t = tiledlayout(numel(rmsAll),numel(corrAll),'TileSpacing','tight');
    % Plot Map
    for i = 1:numel(rmsAll)
        for j = 1:numel(corrAll)
            % Plot each map
            nexttile(t)
            switch MapType
                case {'Map', 'SpecuMap'}
                    try
                        img = Table{i,j}{1};
                        imagesc(img);
                    end
                case 'Stat'
                    try
                        recorded = LoadRfData(parameters.probe, simuDir1);
                        [~, reconstruction] = GenerateParamRecon(recorded);
                        plot(reconstruction.Xmm, Table.linearROI{i,j}{1});
                        ylim([0, 1]);
                    end
                otherwise
                     warning('Unexpected Map type. No data ploted.')
            end
            % Remove x and y axes
            set(gca, 'xtick', [], 'ytick', []);
            end
    end
    % Set axis
    if i == 1
        title(num2str(corrAll(j)), 'FontSize', 8, 'FontWeight', 'normal');
    end
    if j == 1
        ylabel(num2str(corrAll(i)), 'FontSize', 8);
    end
    title(t, 'Correlation length (mm)', 'FontSize', 10)
    ylabel(t, 'Root mean square (mm)', 'FontSize', 10)
    xlabel(t, MapTitle, 'FontSize', 15, 'FontWeight', 'bold')
end 

 