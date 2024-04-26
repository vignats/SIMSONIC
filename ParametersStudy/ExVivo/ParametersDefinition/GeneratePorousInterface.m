bone = '245D';
num_slice = 800;

dirname = ['/calculSSD//Dossier partagé image os exvivo/', bone, '/'];

file = ['SAMPLE_', bone, '_SLICE_', sprintf('%04d', num_slice), '.bmp']; 
filename = [dirname, file];

segmented = false;
boundaryEndost = ExtractBoundary(filename, segmented);
% BOUNDARIES COMPUTATION

bone_bmp = imread(filename); 
threshold = graythresh(bone_bmp); % Find an automatic threshold
binaryImage = imbinarize(bone_bmp, threshold);

boundaries = bwboundaries(binaryImage, 8);
% boundaryEndost = flip(boundaries{1}, 2); %ExtractPorousBoundary(filename, segmented);

figure
imshow(binaryImage);

for k = 1:length(boundaries)
   hold on;
   boundary = boundaries{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end
xlabel('Width');
ylabel('Depth');
grid on
title('Boundary of the endostium'); hold off;


for k = 2:length(boundaries)
    boundary = flip(boundaries{k}, 2);
    if boundary()

function [boundaryEndost] = ExtractPorousBoundary(filename, segmented)
%   This function extracts the boundary of the endost from a binary image.
%
%   Input arguments:
%   - filename: ;bmp image that contain the bone structure.
%
%   Output arguments:
%   - boundaryEndost: Coordinates of the extracted boundary points.

    bone_bmp = imread(filename); 
    
    % Threshold to convert to binary image, if segemented is indicated, the
    % image is already binarize
    if ~segmented
        threshold = graythresh(bone_bmp); % Find an automatic threshold
        binaryImage = imbinarize(bone_bmp, threshold);
    end

    % Select boundary of the endost and periost    
    boundaries = bwboundaries(binaryImage, 8);
    boundary = flip(boundaries{1}, 2);      % First column X-axis & second column Z-axis
    
    % BOUNDARIES COMPUTATION
    % Delimitation of the endost
    % Deletion of the periost by selecting the maximum value along each column.
    [X_bound, idx] = unique(boundary(:,1), 'last');        % Get boundary index and X value for each column. 
    Z_bound = boundary(idx, 2);          % Get corresponfing Z value boundary(idx, 2) and create a matrix of the endost coordinate

    % Deletion of the image boundary
    [nbOccurences, Zcoordinate] = groupcounts(Z_bound);                 % Count the number of occurences to detect the image boundary
    ImageBoundary = Zcoordinate(nbOccurences == max(nbOccurences));     % Define the Z coordinate of the image boundary
    
    Limite = (X_bound(Z_bound == ImageBoundary));                     % Define and center the X coordinate of the image boundary
    LimEndoMax = min(Limite(Limite - mean(Limite) > 0));   
    LimEndoMin = max(Limite(Limite - mean(Limite) < 0));

    if isempty(LimEndoMax)
        LimEndoMax = max(X_bound);
    end

    boundaryEndost = [X_bound(LimEndoMin : LimEndoMax), Z_bound(LimEndoMin : LimEndoMax)];
    boundaryEndost = boundaryEndost';
end

function[] = GetPores(boundaries, boundaryEndost)





end