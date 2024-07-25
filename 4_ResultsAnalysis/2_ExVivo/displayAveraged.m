% This code aim to plot the specular map and the averaged map
clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% LOAD DATAS
dirPath = '~/Documents/BoneRugosity/4_ResultsAnalysis/2_ExVivo/Datas';
bones = {'245D', '227G', '267G'};
zones = {'Z1', 'Z2', 'Z3', 'Z4'};
load(fullfile(dirPath, 'boneSpeed.mat'), 'boneSpeed')
orientation = {{2, 1, 1, 1}, {1, 1, 2, 3}, {1, 2, 1, 1}};

% b = 1; z = 1;
b = 2; z = 3;
% b = 3; z = 1;

%% Get bone averaged image
bone_bmp = imread(fullfile(dirPath, bones{b}, sprintf('AVG_%s_%s_ajustee.bmp', bones{b}, zones{z}))); 

%% Get datas specular
% Get specular map
pathFig = fullfile(dirPath, bones{b}, sprintf('specularity_%s_%s.fig', bones{b}, zones{z}));
fig = openfig(pathFig, 'invisible');

axesObjs = findobj(fig, 'Type', 'axes');       
surfaceObj = findobj(axesObjs(1), 'Type', 'surface');
specularMap= {get(surfaceObj, 'XData'), ...
    get(surfaceObj, 'YData'), get(surfaceObj, 'CData')};

close(fig);

% Get beamformed image
pathFig = fullfile(dirPath, bones{b}, sprintf('Beamform_%s_%s.fig', bones{b}, zones{z}));
fig = openfig(pathFig, 'invisible');

axesObjs = findobj(fig, 'Type', 'axes');
childrenObjs = get(axesObjs(1), 'Children');
beamformMap{1} = {get(childrenObjs, 'XData'), ...
    get(childrenObjs, 'YData'), get(childrenObjs, 'CData')};

close(fig);

% Get endosteum segmentation
load(fullfile(dirPath, bones{b}, sprintf('Parabola_%s_%s.mat', bones{b}, zones{z})));
rep = orientation{b}{z};
segmentedEndost{1} = [X_ENDO{rep}; Z_ENDO{rep}].*1000;
segmentedPeriost{1} = [X_PERI{rep}; Z_PERI{rep}].*1000;

%% Get boundaries
speedSound = getfield(boneSpeed, sprintf('Bone%s', bones{b}), zones{z});
boundaryWidth = speedSound{1}/(2*2.6e3);
[surface, surfaceMax] = ExpandParabola(segmentedEndost{1}, boundaryWidth, 'specularity');
[~, surfaceMin] = ExpandParabola(segmentedEndost{1}, -boundaryWidth, 'specularity');

%% Plot figures
stepXRay = 0.009; % Step of the X-Ray image (mm)
X = 0:stepXRay:size(bone_bmp, 2)*stepXRay; X = X - mean(X);
Z = flipud(0:stepXRay:size(bone_bmp, 1)*stepXRay);

figure
ax(1) = subplot(1,2,1);
imagesc(X, Z, bone_bmp)
xlabel('Width (mm)', 'Interpreter', 'latex',  'FontSize', 16);
ylabel('Depth (mm)', 'Interpreter', 'latex',  'FontSize', 16);
title('Ex-vivo Bone Averaged Image', 'Interpreter', 'latex',  'FontSize', 18)
axis image
colormap(ax(1),gray)
ax(2) = subplot(1,2,2);
imagesc(specularMap{1}, specularMap{2}-2, specularMap{3})
hold on 
plot(surfaceMax(1,:), surfaceMax(2,:)-2, 'r', 'LineWidth', 1)
plot(surfaceMin(1,:), surfaceMin(2,:)-2, 'r', 'LineWidth', 1)
colormap(ax(2),jet)
colorbar
xlabel('Width (mm)', 'Interpreter', 'latex',  'FontSize', 16);
ylabel('Depth (mm)', 'Interpreter', 'latex',  'FontSize', 16);
title('Specularity Map of the Bone', 'Interpreter', 'latex',  'FontSize', 18)

axis image
grid off
clim([0 1])
