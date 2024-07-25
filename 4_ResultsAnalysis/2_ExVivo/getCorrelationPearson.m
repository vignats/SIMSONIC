% This code aim to post-process the ex-vivo specualrity maps
clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% Load datas
dirPath = '~/Documents/BoneRugosity/4_ResultsAnalysis/2_ExVivo/Datas';
load(fullfile(dirPath, 'experimentalData.mat'), 'meanNum', 'meanRoi');
% The results can be boneParameters or boneParametersKc
load(fullfile(dirPath, 'boneParametersKc.mat'));

%% To delete '271G'
fields = fieldnames(boneParametersKc);
for i = 1:numel(fields)
    field = fields{i};
    % Delete the last row corresponding to 271G
    boneParametersKc(1).(field) = boneParametersKc(1).(field)(1:end-1, :);
end

meanNum = meanNum(1:end-1, :);
meanRoi = meanRoi(1:end-1, :);
%% Correlation between bone parameters and specularity metrics
% Outputs and inputs datas in column
outputRoi = table2array(meanRoi); outputRoi = outputRoi(:);
outputNum = table2array(meanNum); outputNum = outputNum(:);

% The results can be boneParametersKc or boneParametersKc
iRq = table2array(boneParametersKc.Rq); iCorr = table2array(boneParametersKc.Corr);
iEP = table2array(boneParametersKc.EPore); idP = table2array(boneParametersKc.dPore);
inputParameters = [iRq(:) iCorr(:) iEP(:) idP(:)];

[corrRoi, PvaluesR] = corrcoef([inputParameters outputRoi]);
[corrNum, PvaluesN] = corrcoef([inputParameters outputNum]);

%% Plot datas
Xname = {'Rq (mm)', 'Correlation length (mm)', 'Porosity (%)', 'Pore diameter (mm)'};
col = {'r', 'b', 'g', 'y'};
bones = {'245D', '227G', '267G', '271G'};
data = {'Rq', 'Corr', 'EPore', 'dPore'};

figure;
t = tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(fieldnames(boneParametersKc))
    nexttile(t)
    input = table2array(getfield(boneParametersKc, data{i}));
    for j = 1:size(meanRoi,1)
        scatter(input(j,:), table2array(meanRoi(j,:)), '*', 'MarkerEdgeColor', col{j})
        hold on
        title(Xname{i})
    end
    if i == 1
        ylabel('Mean specularity in the Region of interest (%)')
        legend(bones)
    end
    hold off
end

for i = 1:numel(fieldnames(boneParametersKc))
    nexttile(t)
    input = table2array(getfield(boneParametersKc, data{i}));
    for j = 1:size(meanRoi,1)
        scatter(input(j,:), table2array(meanNum(j,:)), '*', 'MarkerEdgeColor', col{j})
        hold on
        title(Xname{i})
    end
    if i == 1
        ylabel('Specular pixel in the Region of interest (%)')
    end
    hold off
end

%% Plot interesting ones
figure;
t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile(t)
i = 2;
input = table2array(getfield(boneParametersKc, data{i}));
for j = 1:size(meanRoi,1)
    scatter(input(j,:), table2array(meanNum(j,:)), 'd', 'filled')
    hold on
    xlabel(Xname{i}, 'Interpreter', 'latex',  'FontSize', 18)
end
ylabel('Specular pixel in the Region of interest (\%)', 'Interpreter', 'latex',  'FontSize', 18)
ax = gca;
set(ax, 'Units', 'normalized'); 
ax.FontSize = 16;

nexttile(t)
i = 4;
input = table2array(getfield(boneParametersKc, data{i}));
for j = 1:size(meanRoi,1)
    scatter(input(j,:), table2array(meanNum(j,:)), 'd', 'filled')
    hold on
    xlabel(Xname{i}, 'Interpreter', 'latex',  'FontSize', 18)
end
legend(bones, 'Interpreter', 'latex',  'FontSize', 16)
ax = gca;
set(ax, 'Units', 'normalized'); 
ax.FontSize = 16;