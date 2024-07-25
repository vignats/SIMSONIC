% COMPUTE PARAMETERS OF EX-VIVO BONES ON ONE SLICES PER ZONE
clear; 
close all;
clc;
addpath(genpath('~/Documents')); 
addpath(genpath('/calculSSD/salome')); 

%% Load parameters
bones = {'245D', '227G', '267G', '271G'};
slices = {{1134, 3195, 3852, 5511}, {2002, 3595, 5614, 6721}, {1700, 3410, 5478, 6716}}; 
kcAll = [0.041 0.046 0.057 0.064; 0.051 0.046 0.033 0.037; 0.033 0.027 0.030 0.027];
load('/calculSSD/salome/Simulation-10mai/boneSpeed.mat', 'boneSpeed')

%% Display figures
figure
t = tiledlayout(numel(bones), numel(slices{1}), 'Padding', 'compact', 'TileSpacing', 'compact');

for b = 1:numel(bones)  
    bone.id = bones{b};         % Bone from ex-vivo files
    for s = 1:numel(slices{b})
        bone.image = slices{b}{s};        % Slice selected

        % Load data 
        simulation_name = sprintf('Bone%s-Image%04d/', bone.id, bone.image);
        load(fullfile('/calculSSD/salome/Simulation-10mai/', simulation_name, 'parameters.mat'));

        % Binarize bone image
        dirname = '~/Documents/BoneRugosity/2_ParametersDefinition/BoneImage';
        file = sprintf('SAMPLE_%s_SLICE_%04d.bmp', bone.id, bone.image); 
        filename = fullfile(dirname, file); 

        bone_bmp = imread(filename); 
        threshold = graythresh(bone_bmp); % Find an automatic threshold
        binaryImage = imbinarize(bone_bmp, threshold);

        if b == 1 && s ~=1
            if s == 2 || s == 4
                angle = -10;
            elseif s == 3
                angle = -15;
            end
            binaryImage = imrotate(binaryImage, angle,'bilinear','loose');
            binaryImage = binaryImage(1:find(...
            sum(binaryImage(: , 1: size(binaryImage, 2)/2),2) > 10, 1, 'last'), :);
        end 

        % Extract endost
        simuName = fullfile('/calculSSD/salome/Simulation-10mai/', simulation_name);
        lambda =  boneSpeed.(simuName(36:43)).(simuName(45:53)){1} / (probe.fc)*1e3; %(medium.cp(2) / (probe.fc))*1e3; 
        [profile, roughness, ~, kc] = FilterProfil(binaryImage, kcAll(b,s));
        
        % Compute the correlation length and rms of the roughness profile 
        Rms = rms(roughness);
        corr = ComputeCorr(roughness, grid.step);

        % Compute the porosity in a surface of one wavelength around the
        % boundary
        speedSound = getfield(boneSpeed, sprintf('Bone%s', bones{b}), sprintf('Image%04d', slices{b}{s}));
        boundaryWidth = speedSound{1}/(2.6e3);    %Frequency of the probe : 2.6MHz
        
        [boundaryEndost, boundaryPores] = ExtractBoundary(binaryImage); % Extract the pores
        [porosity, poreSize, imageDisplay] = ComputePorosity(binaryImage, boundaryEndost, boundaryPores, boundaryWidth, false);

        % Plot results
        nexttile(t)
        stepXRay = 0.009; % Step of the X-Ray image (mm)
        X = 0:stepXRay:size(imageDisplay, 2)*stepXRay; %X = X - mean(X);
        Z = flipud(0:stepXRay:size(imageDisplay, 1)*stepXRay);
        imagesc(X, Z, imageDisplay)
        hold on
        plot(boundaryEndost(1,:)*stepXRay, boundaryEndost(2,:)*stepXRay, 'r', 'LineWidth', 2)
        xlabel('Width (mm)');
        ylabel('Depth (mm)');
        hold off
        axis image

        % Calculate dimensions of current subplot
        ax = gca;
        set(ax, 'Units', 'normalized'); 
        pos = ax.Position;

        % Adjust position of annotation box
        dim = [pos(1), pos(2) + pos(4) - 0.08, 0.3, 0.1];
        str = {[sprintf('Rq = %.2fmm', Rms)]...
            [texlabel('rho') sprintf(' = %.3fmm', corr)] sprintf('E.Por = %.1f%%', porosity) sprintf('d.Por = %.3fmm', poreSize)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on', 'BackgroundColor','w',...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    end
end
title(t, 'Parameters of Bone Health on Ex-Vivo X-Ray Images', 'Interpreter', 'latex',  'FontSize', 22)
