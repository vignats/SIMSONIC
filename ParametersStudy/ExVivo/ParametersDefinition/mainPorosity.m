clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% Create bone structure 
% bones = {'271G' '251G' '245D' '227G'};
% min_slices = [1590, 761, 470, 1210];
% max_slices = [2954, 2125, 1861, 2793];
% % boneTable = table(min_slices', max_slices', 'VariableNames', {'min_slice', 'max_slice'}, 'RowNames', bones);

%% Get binarize bone image
bone = '245D';
num_slice = 800;

dirname = ['/calculSSD//Dossier partagé image os exvivo/', bone, '/'];
file = ['SAMPLE_', bone, '_SLICE_', sprintf('%04d', num_slice), '.bmp']; 
filename = fullfile(dirname, file);

bone_bmp = imread(filename); 
threshold = graythresh(bone_bmp); % Find an automatic threshold
binaryImage = imbinarize(bone_bmp, threshold);

%% Compute the porosity of the bone
% Compute X-Ray image and ultrasound rf informations 
[boundaryEndost, boundaryPores] = ExtractBoundary(filename, false);

% Compute the porosity in a surface of one wavelength around the boundary
nbWavelength = 1/2;
porosity = ComputePorosity(binaryImage, boundaryEndost, nbWavelength, true);

