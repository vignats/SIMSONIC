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
bone = '267G';
zone = 2;
num_slice = 3410;

% dirname = '/calculSSD/salome/BoneImage/';
dirname = ['/calculSSD/Dossier partagé image os exvivo/', bone, '/ZONE_US_0', int2str(zone)];

file = ['SAMPLE_', bone, '_SLICE_', sprintf('%08d', num_slice), '.bmp'];

filename = fullfile(dirname, file);

bone_bmp = imread(filename); 
threshold = graythresh(bone_bmp); % Find an automatic threshold
binaryImage = imbinarize(bone_bmp, threshold);

%% Compute the porosity of the bone
% Compute X-Ray image and ultrasound rf informations 
[boundaryEndost, boundaryPores] = ExtractBoundary(binaryImage);

% Compute the porosity in a surface of one wavelength around the boundary
width = 1/2;
porosity = ComputePorosity(binaryImage, boundaryEndost, width, true);

