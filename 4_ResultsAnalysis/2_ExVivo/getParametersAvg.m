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

boneParametersAvg = struct('Rq', type, 'Corr', type, 'EPore', type,...
    'dPore', type, 'kc', type);

for b = 1:numel(bones)  
    for z = 1:numel(zones)        
        % LOAD DATAS
        bone_bmp = imread(fullfile(dirPath, bones{b}, sprintf('AVG_%s_Z%01d.bmp', bones{b}, z)));
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
        [~, roughness, ~, boneParametersAvg.kc{b,z}] = FilterProfil(binaryImage,0);
        
        % Compute the correlation length and rms of the roughness profile 
        boneParametersAvg.Rq{b,z} = rms(roughness);
        boneParametersAvg.Corr{b,z} = ComputeCorr(roughness, 0.009);    % X-Ray image resolution : 9um

        % Compute the porosity in a surface of one wavelength above the boundary
        speedSound = getfield(boneSpeed, sprintf('Bone%s', bones{b}), zones{z});
        boundaryWidth = speedSound{1}/(2.6e3);    %Frequency of the probe : 2.6MHz
        
        [boundaryEndost, boundaryPores] = ExtractBoundary(binaryImage); % Extract the pores
        [porosity, poreSize, displayImage] = ComputePorosity(binaryImage, boundaryEndost, boundaryPores, boundaryWidth, false);
            
        boneParametersAvg.EPore{b,z} = porosity;
        boneParametersAvg.dPore{b,z} = poreSize;
    end
end
        