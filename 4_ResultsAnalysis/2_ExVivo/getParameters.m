% This code aim to compute the parameters over each zone of each bone
% samples
clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% Load datas
dirPath = '~/Documents/BoneRugosity/4_ResultsAnalysis/2_ExVivo/Datas';
load(fullfile(dirPath, 'boneSpeed.mat'), 'boneSpeed')

pathFigAll = '/calculSSD/Dossier partagé image os exvivo/';
bones = {'245D', '227G', '267G', '271G'};
zones = {'Z1', 'Z2', 'Z3', 'Z4'};
kcAll = {{0.041, 0.046, 0.057, 0.064}, {0.051, 0.046, 0.033, 0.037},...
    {0.033, 0.027, 0.030, 0.027}, {0, 0, 0, 0}};

%% Compute bone microstructure parameters
type = table('Size', [numel(bones), numel(zones)], ...
                'VariableType', repmat({'double'}, 1, numel(zones)), ...
                'VariableNames', cellstr(string(zones)), ...
                'RowNames', cellstr(string(bones)));

meanParameters = struct('Rq', type, 'Corr', type, 'EPore', type,...
    'dPore', type, 'kc', type);

for b = 1:numel(bones)  
    for z = 1:numel(zones)        
        % LOAD DATAS
        pathFigure = fullfile(pathFigAll, bones{b}, sprintf('ZONE_US_%02d', z));
        figAll = dir(fullfile(pathFigure, 'SAMPLE*.bmp'));
        
        % COMPUTE PARAMETERS 
        Rq = 0; Corr = 0; EPore = 0; dPore = 0; numFig = 0; kc = 0;
        for idx = 1:10:numel(figAll)
            try
                % Binarize bone image
                bone_bmp = imread(fullfile(pathFigure, figAll(idx).name)); 
                threshold = graythresh(bone_bmp); % Find an automatic threshold
                binaryImage = imbinarize(bone_bmp, threshold);
                
                % Rotate image if necessary
                if b == 1 && z ~=1
                    if z == 2 || z == 4
                        angle = -10;
                    elseif z == 3
                        angle = -15;
                    end
                    binaryImage = imrotate(binaryImage, angle,'bilinear','loose');
                    binaryImage = binaryImage(1:find(...
                    sum(binaryImage(: , 1: size(binaryImage, 2)/2),2) > 10, 1, 'last'), :);
                end
    
                % Extract endosteum
                [profile, roughness, xP, kcSlice] = FilterProfil(binaryImage,0);
                
                % Compute the correlation length and rms of the roughness profile 
                Rq = Rq + rms(roughness);
                Corr = Corr + ComputeCorr(roughness, 0.009);    % X-Ray image resolution : 9um
                kc = kc + kcSlice;

                % Compute the porosity in a surface of one wavelength above the boundary
                speedSound = getfield(boneSpeed, sprintf('Bone%s', bones{b}), zones{z});
                boundaryWidth = speedSound{1}/(2.6e3);    %Frequency of the probe : 2.6MHz
                
                [boundaryEndost, boundaryPores] = ExtractBoundary(binaryImage); % Extract the pores
                [porosity, poreSize, displayImage] = ComputePorosity(binaryImage, boundaryEndost, boundaryPores, boundaryWidth, false);
                EPore = EPore + porosity;
                dPore = dPore + poreSize;
    
                numFig = numFig + 1;
            catch
                fprintf('Did not worked for bone %s in %s for the image %s \n', bones{b}, zones{z}, figAll(idx).name) 
            end
        end
        
        % Compute the mean parameters
        meanParameters.Rq{b,z} = Rq / numFig;
        meanParameters.Corr{b,z} = Corr / numFig;
        meanParameters.EPore{b,z} = EPore / numFig;
        meanParameters.dPore{b,z} = dPore / numFig;
        meanParameters.kc{b,z} = kc/numFig;

        fprintf('Number of samples for bone %s = %d \n', bones{b},numFog)
    end
end

