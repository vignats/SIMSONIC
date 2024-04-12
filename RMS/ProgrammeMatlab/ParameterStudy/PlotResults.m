clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% PROCESSING OF THE RF DATA 
% Generate a table to stock the geometry and the specularity Map
pathSimuAll = '/calculSSD/salome/Simulation-04avr';
rmsAll = 0.03 + (0:9) * 0.05;      % List of rms values (mm)
corrAll = [0.5 1 2 4];             % List of correlation length (mm)

SpecularAll = struct();

%% COMPUTE GEOMETRY MAP
SpecularAll.Map = ComputeMap('Map', rmsAll, corrAll, pathSimuAll, SpecularAll);

%% PLOT GEOMETRY MAPS
MapTitle = 'Interface profile for various RMS and correlation length';
PlotAll('Map', SpecularAll, rmsAll, corrAll, pathSimuAll, MapTitle);

%% COMPUTE SPECULARITY MAPS 
SpecularAll.SpecuMap = ComputeMap('SpecuMap', rmsAll, corrAll, pathSimuAll, SpecularAll);

%% PLOT SPECULARITY MAPS
MapTitle = 'Specular probability for various RMS and correlation length';
PlotAll('SpecuMap', SpecularAll, rmsAll, corrAll, pathSimuAll, MapTitle);

%% COMPUTE CORRESPONDING STATISTICAL INFORMATIONS
SpecularAll.SpecuProba = ComputeMap('SpecuProba', rmsAll, corrAll, pathSimuAll, SpecularAll);

%% PLOT ALL MAPS
MapTitle = 'Specular probability along the lateral position';
PlotAll('SpecuProba', SpecularAll, rmsAll, corrAll, pathSimuAll, MapTitle);

%% TO USE IN EXCEL 
[meanROI, stdROI, corrROI] = deal(table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'double'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll))));

for i = 1:numel(rmsAll)
    for j = 1:numel(corrAll)
        try
            meanROI{i,j} = SpecularAll.SpecuProba{i,j}{1}.meanROI;
            stdROI{i,j} = SpecularAll.SpecuProba{i,j}{1}.stdROI;
            % corrROI{i,j} = SpecuProbaAll{i,j}{1}.corrROI;
        end
    end
end

%% TO SUBPLOT
format = 'simulation_rms_%.2f_cl_%.1f/';
simuDir1 = fullfile(pathSimuAll, sprintf(format, rmsAll(1), corrAll(1)));
parameters = load(fullfile(simuDir1, 'parameters.mat'));
recorded = LoadRfData(parameters.probe, simuDir1);
[~, reconstruction] = GenerateParamRecon(recorded);

legendPlot = {};
figure
for i = 1:numel(rmsAll)
    for j = 1:numel(corrAll)
        if ~isempty(SpecularAll.SpecuProba{i,j}{1})
            plot(reconstruction.Xmm, SpecularAll.SpecuProba{i,j}{1}.linearROI);
            ylim([0, 1]);
            legendPlot{end+1} = ['RMS = ', SpecularAll.SpecuProba.Properties.RowNames{i}, ' CORR = ', SpecularAll.SpecuProba.Properties.VariableNames{j}];

            hold on
        end
    end
end
xlabel('Lateral position (mm)', Interpreter='latex')
ylabel('Specular probability', Interpreter='latex')
title('Specular probability along the lateral position');
legend(legendPlot)
ylim([0, 1]);

function[Table] = ComputeMap(MapType, rmsAll, corrAll, pathSimuAll, SpecularAll)
    % Compute parameters once 
    format = 'simulation_rms_%.2f_cl_%.1f/';
    simuDir1 = fullfile(pathSimuAll, sprintf(format, rmsAll(1), corrAll(1)));
    parameters = load(fullfile(simuDir1, 'parameters.mat'));
    recorded = LoadRfData(parameters.probe, simuDir1);
    [~, reconstruction] = GenerateParamRecon(recorded);

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
                case 'SpecuProba'
                    try
                        probability = ProbaROI(SpecularAll.SpecuMap{i,j}{1}, reconstruction, parameters, false);
                        Table{i,j}{1} = probability;
                        disp(i)
                    end
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
    recorded = LoadRfData(parameters.probe, simuDir1);
    [~, reconstruction] = GenerateParamRecon(recorded);

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
                case 'SpecuProba'
                    try
                        plot(reconstruction.Xmm, Table{i,j}{1}.linearROI);
                        ylim([0, 1]);
                    end
                otherwise
                     warning('Unexpected Map type. No data ploted.')
            end
            % Remove x and y axes
            set(gca, 'xtick', [], 'ytick', []);
            
            % Set axis
            if i == 1
                title(Table.Properties.VariableNames{j}, 'FontSize', 8, 'FontWeight', 'normal');
            end
            if j == 1
                ylabel(Table.Properties.RowNames{i}, 'FontSize', 8);
            end
        end
    end
    title(t, 'Correlation length (mm)', 'FontSize', 10)
    ylabel(t, 'Root mean square (mm)', 'FontSize', 10)
    xlabel(t, MapTitle, 'FontSize', 15, 'FontWeight', 'bold')
end 

 