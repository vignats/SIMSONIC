% This code aim to post-process the experimental ex-vivo specualrity maps
clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% LOAD DATAS
dirPath = '~/Documents/BoneRugosity/4_ResultsAnalysis/2_ExVivo/Datas';
bones = {'245D', '227G', '267G', '271G'};
zones = {'Z1', 'Z2', 'Z3', 'Z4'};
load(fullfile(dirPath, 'boneSpeed.mat'), 'boneSpeed')
orientation = {{2, 1, 1, 1}, {1, 1, 2, 3}, {1, 2, 1, 1}};

%% Create
specularMap = table('Size', [numel(bones), numel(zones)], ...
                'VariableType', repmat({'cell'}, 1, numel(zones)), ...
                'VariableNames', cellstr(string(zones)), ...
                'RowNames', cellstr(string(bones)));

beamformMap = table('Size', [numel(bones), numel(zones)], ...
                'VariableType', repmat({'cell'}, 1, numel(zones)), ...
                'VariableNames', cellstr(string(zones)), ...
                'RowNames', cellstr(string(bones)));

segmentedEndost = table('Size', [numel(bones), numel(zones)], ...
                'VariableType', repmat({'cell'}, 1, numel(zones)), ...
                'VariableNames', cellstr(string(zones)), ...
                'RowNames', cellstr(string(bones)));

segmentedPeriost = table('Size', [numel(bones), numel(zones)], ...
                'VariableType', repmat({'cell'}, 1, numel(zones)), ...
                'VariableNames', cellstr(string(zones)), ...
                'RowNames', cellstr(string(bones)));

%% COMPLETE TABLES
for b = 1:length(bones)
    for z = 1:length(zones)
        % Get specularity map
        pathFig = fullfile(dirPath, bones{b}, sprintf('specularity_%s_%s.fig', bones{b}, zones{z}));
        fig = openfig(pathFig, 'invisible');

        axesObjs = findobj(fig, 'Type', 'axes');       
        surfaceObj = findobj(axesObjs(1), 'Type', 'surface');
        specularMap{b,z}{1} = {get(surfaceObj, 'XData'), ...
            get(surfaceObj, 'YData'), get(surfaceObj, 'CData')};

        close(fig);

        % Get beamformed image
        pathFig = fullfile(dirPath, bones{b}, sprintf('Beamform_%s_%s.fig', bones{b}, zones{z}));
        fig = openfig(pathFig, 'invisible');

        axesObjs = findobj(fig, 'Type', 'axes');
        childrenObjs = get(axesObjs(1), 'Children');
        if length(childrenObjs) == 1 % If the segmented endosteum isn't in the figure
            beamformMap{b,z}{1} = {get(childrenObjs, 'XData'), ...
            get(childrenObjs, 'YData'), get(childrenObjs, 'CData')};
        
            close(fig);
    
            % Get endosteum segmentation
            load(fullfile(dirPath, bones{b}, sprintf('Parabola_%s_%s.mat', bones{b}, zones{z})));
            rep = orientation{b}{z};
            segmentedEndost{b,z}{1} = [X_ENDO{rep}; Z_ENDO{rep}].*1000;
            segmentedPeriost{b,z}{1} = [X_PERI{rep}; Z_PERI{rep}].*1000;
        else
            beamformMap{b,z}{1} = {get(childrenObjs(3), 'XData'), ...
            get(childrenObjs(3), 'YData'), get(childrenObjs(3), 'CData')};
    
            % Get endosteum segmentation
            segmentedEndost{b,z}{1} = [get(childrenObjs(1), 'XData'); ...
            get(childrenObjs(1), 'YData')];
            segmentedPeriost{b,z}{1} = [get(childrenObjs(2), 'XData'); ...
            get(childrenObjs(2), 'YData')];
            close(fig);
        end
        
    end
end

%% PLOT SPECULAR IMAGES
figure
t = tiledlayout(numel(bones), numel(zones), 'Padding', 'compact', 'TileSpacing', 'compact');
for b = 1:length(bones)
    for z = 1:length(zones)
        % Extract boundaries
        speedSound = getfield(boneSpeed, sprintf('Bone%s', bones{b}), zones{z});
        boundaryWidth = speedSound{1}/(2*2.6e3);
        [surface, surfaceMax] = ExpandParabola(segmentedEndost{b,z}{1}, boundaryWidth, 'specularity');
        [~, surfaceMin] = ExpandParabola(segmentedEndost{b,z}{1}, -boundaryWidth, 'specularity');
        
        % Plot figures
        % figure
        nexttile(t)
        imagesc(specularMap{b,z}{1}{1}, specularMap{b,z}{1}{2}, specularMap{b,z}{1}{3})
        hold on 
        % plot(segmentedEndost{b,s}{1}(1,:), segmentedEndost{b,s}{1}(2,:), 'r', 'LineWidth', 3)
        % plot(segmentedPeriost{b,s}{1}(1,:), segmentedPeriost{b,s}{1}(2,:), 'b', 'LineWidth', 3)
        plot(surfaceMax(1,:), surfaceMax(2,:), 'r', 'LineWidth', 3)
        plot(surfaceMin(1,:), surfaceMin(2,:), 'r', 'LineWidth', 3)
        % plot(surface(:,1), surface(:,2), 'w', 'LineWidth', 3)
        hold off
        colormap jet
        colorbar
        axis image
        grid off
        clim([0 1])
    end
end

%% COMPUTE RESULTS
meanRoi = table('Size', [numel(bones), numel(zones)], ...
                'VariableType', repmat({'double'}, 1, numel(zones)), ...
                'VariableNames', cellstr(string(zones)), ...
                'RowNames', cellstr(string(bones)));

meanNum = table('Size', [numel(bones), numel(zones)], ...
                'VariableType', repmat({'double'}, 1, numel(zones)), ...
                'VariableNames', cellstr(string(zones)), ...
                'RowNames', cellstr(string(bones)));

meanSpecu = table('Size', [numel(bones), numel(zones)], ...
                'VariableType', repmat({'double'}, 1, numel(zones)), ...
                'VariableNames', cellstr(string(zones)), ...
                'RowNames', cellstr(string(bones)));

for b = 1:length(bones)
    for z = 1:length(zones)
        % Create tables to store results 
        meanValues = zeros(1, length(surfaceMax));
        sumValues = zeros(2, length(surfaceMax));
        meanS = zeros(2, length(surfaceMax));

        % Extract boundaries
        speedSound = getfield(boneSpeed, sprintf('Bone%s', bones{b}), zones{z});
        boundaryWidth = speedSound{1}/(2*2.6e3);
        [~, surfaceMax] = ExpandParabola(segmentedEndost{b,z}{1}, boundaryWidth, 'specularity');
        [~, surfaceMin] = ExpandParabola(segmentedEndost{b,z}{1}, -boundaryWidth, 'specularity');
        
        % Compute mean specular index 
        for i = 1:length(surfaceMax)
            % Extract index of the ROI
            indSpecu = find(specularMap{b,z}{1}{1} >= surfaceMax(1,i), 1);
            indMax = find(specularMap{b,z}{1}{2} >= surfaceMax(2,i), 1);
            indMin = find(specularMap{b,z}{1}{2} <= surfaceMin(2,i), 1, 'last');

            meanValues(i) = mean(specularMap{b,z}{1}{3}(indMax:indMin, indSpecu));
            meanS(i) = mean(specularMap{b,z}{1}{3}(indMax:indMin, indSpecu) >= 0.4);
            sumValues(1, i) = sum(specularMap{b,z}{1}{3}(indMax:indMin, indSpecu) >= 0.4);
            sumValues(2, i) = numel(specularMap{b,z}{1}{3}(indMax:indMin, indSpecu));
        end

        meanRoi{b,z} = mean(meanValues).*100;
        meanSpecu{b,z} = mean(meanS, 'all').*100;
        meanNum{b,z} = sum(sumValues(1,:))/sum(sumValues(2,:)).*100;
    end
end

%% COMPUTE HISTOGRAMME
figure
t = tiledlayout(numel(bones), numel(zones), 'Padding', 'compact', 'TileSpacing', 'compact');
for b = 1:length(bones)
    for z = 1:length(zones)
         % Extract boundaries
        speedSound = getfield(boneSpeed, sprintf('Bone%s', bones{b}), zones{z});
        boundaryWidth = speedSound{1}/(2*2.6e3);
        [~, surfaceMax] = ExpandParabola(segmentedEndost{b,z}{1}, boundaryWidth, 'specularity');
        [~, surfaceMin] = ExpandParabola(segmentedEndost{b,z}{1}, -boundaryWidth, 'specularity');
        
        pixelValues = [];
        for i = 1:length(surfaceMax)
            % Compute index of the ROI position
            indSpecu = find(specularMap{b,z}{1}{1} >= surfaceMax(1,i), 1);
            indMax = find(specularMap{b,z}{1}{2} >= surfaceMax(2,i), 1);
            indMin = find(specularMap{b,z}{1}{2} <= surfaceMin(2,i), 1, 'last');

            % Values of the pixels
            pixelValues = [pixelValues; specularMap{b,z}{1}{3}(indMax:indMin, indSpecu)];
        end
        nexttile(t)
        histogram(pixelValues, 10)
    end
end
%% VERIFY SEGMENTATION
b = 3;
z = 4;

load(fullfile(dirPath, bones{b}, sprintf('Parabola_%s_%s.mat', bones{b}, zones{z})));
figure
t = tiledlayout(2, 5, 'Padding', 'compact');
for rep = 1:10
    % % Beamformed maps
    % pathFig = fullfile('/calculSSD/Dossier partagé image os exvivo/specularity_from_amadou', ...
    %     bones{b}, ['specular_only_1D_exact_Phantom_267G_B', sprintf('%02d', rep), '_1_wind_length_050mod', '.fig']);
    % fig = openfig(pathFig, 'invisible');
    % 
    % axesObjs = findobj(fig, 'Type', 'axes');
    % childrenObjs = get(axesObjs(1), 'Children');
    % 
    % beamformMap{b,s}{1} = {get(childrenObjs(3), 'XData'), ...
    %     get(childrenObjs(3), 'YData'), get(childrenObjs(3), 'CData')};
    % 
    % close(fig);
    
    % Specular map
    pathFig = fullfile('/calculSSD/Dossier partagé image os exvivo/specularity_from_amadou',...
        bones{b}, sprintf('specularity_1D_exact_Phantom_%s_T%02d_1.fig', bones{b}, rep));
    fig = openfig(pathFig, 'invisible');

    axesObjs = findobj(fig, 'Type', 'axes');       
    surfaceObj = findobj(axesObjs(1), 'Type', 'surface');
    specularMap{b,z}{1} = {get(surfaceObj, 'XData'), ...
        get(surfaceObj, 'YData'), get(surfaceObj, 'CData')};

    close(fig);

    % Plot fig
    nexttile(t)
    imagesc(specularMap{b,z}{1}{1}, specularMap{b,z}{1}{2}, specularMap{b,z}{1}{3})
    hold on 
    plot(X_ENDO{rep}.*1000, Z_ENDO{rep}.*1000, 'r', 'LineWidth', 3)
    plot(X_PERI{rep}.*1000, Z_PERI{rep}.*1000, 'b', 'LineWidth', 3)
    axis image
    colormap jet
    % colormap gray
    % clim([-40 0])
end