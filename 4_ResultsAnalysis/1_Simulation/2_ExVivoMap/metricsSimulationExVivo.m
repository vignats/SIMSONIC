% This code aim to post-process the ex-vivo specualrity maps
clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% Initiating
pathSimuAll = '/calculSSD/salome/Simulation-10mai'; 

bones = {'245D', '227G', '267G'};
slices = {'Z1', 'Z2', 'Z3', 'Z4', 'Z0'};
slicesNb = {{1134, 3195, 3852, 5511, 0000}, {2002, 3595, 5614, 6721, 0000}, {1700, 3410, 5478, 6716, 0000}}; 

meanROI = table('Size', [numel(bones), numel(slices)], ...
                'VariableType', repmat({'double'}, 1, numel(slices)), ...
                'VariableNames', cellstr(string(slices)), ...
                'RowNames', cellstr(string(bones)));

%% Compute the parameters for all the bones slices
for b = 1:length(bones)
    [lowBoundary, highBoundary, X_bound_low, X_bound_high, meanValue] = ExtractROI(bones{b});
    meanROI{b,5} = meanValue;
    for s = 1:length(slices)-1
        postProcess = load(fullfile(pathSimuAll, sprintf('Bone%s-Image%04d/', bone.id, slicesNb{b}{s}), 'postProcess.mat'));
        imageBone = postProcess.SpecularProbaMap.Bone;
        % Mean value inside the ROI
        valueRoi = zeros(min(length(lowBoundary), length(highBoundary)), 1);
        for i = 1: min(length(lowBoundary), length(highBoundary))
            bound = X_bound_low;
            if numel(X_bound_high) == min(length(lowBoundary), length(highBoundary))
                bound = X_bound_high;
            end
            valueRoi = mean(imageBone(lowBoundary(i) : highBoundary(i), bound(i)));
        end
        
        meanROI{b,s} = mean(valueRoi);

        figure
        imagesc(imageBone)
        colormap jet
        hold on
        plot(X_bound_low, lowBoundary, 'g', 'LineWidth', 2)
        hold on
        plot(X_bound_high, highBoundary, 'r', 'LineWidth', 2)
        axis image
    end 
end
%%
[lowBoundary, highBoundary, X_bound_low, X_bound_high, meanValue] = ExtractROI(bones{b});
    meanROI{b,5} = meanValue;
%% Functions
function[lowBoundary, highBoundary, X_bound_low, X_bound_high, meanValue] = ExtractROI(bone)
    pathSimuAll = '/calculSSD/salome/Simulation-10mai';
    postProcess = load(fullfile(pathSimuAll, ['Bone' bone '-Image0000/'], 'postProcess.mat'));
    image = postProcess.SpecularProbaMap.Bone;
    
    % Creation of a mask 
    simu0 = image;
    simu0(simu0 > 0.4) = 1;
    simu0(simu0 < 0.4) = 0;
    
    boundaries = bwboundaries(simu0);
    maxNumElements = 0;
    
    for k = 1: length(boundaries)
        if numel(boundaries{k}) > maxNumElements
            maxNumElements = numel(boundaries{k});
            boundary = boundaries{k};
        end
    end
    
    [~, idx] = unique(boundary(:,2), 'last'); 
    X_bound_high = boundary(min(idx) : max(idx), 2);  
    X_bound_low = boundary(1 : min(idx) - 1, 2);  
    
    highBoundary = boundary(min(idx) : max(idx), 1);    % As 1 pixel = 0.15 mm and lambda = 1.4 mm
    lowBoundary = boundary(1:min(idx) -1, 1);
    
    valueRoi = zeros(min(length(lowBoundary), length(highBoundary)), 1);
    for i = 1: min(length(lowBoundary), length(highBoundary))
        bound = X_bound_low;
        if length(X_bound_high) == min(length(lowBoundary), length(highBoundary))
            bound = X_bound_high;
        end
        valueRoi = mean(image(lowBoundary(i) : highBoundary(i), bound(i)));
    end
    meanValue = mean(valueRoi);
  
    figure
    imagesc(simu0)
    hold on
    % plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
    hold on
    plot(X_bound_high, highBoundary, 'g', 'LineWidth', 2)
    hold on
    plot(X_bound_low, lowBoundary, 'r', 'LineWidth', 2)
    axis image
end