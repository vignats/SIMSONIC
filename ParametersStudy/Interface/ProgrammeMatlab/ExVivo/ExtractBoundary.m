function[boundaryEndost, boundaryPores] = ExtractBoundary(filename, segmented)
%   This function extracts the boundary of the endost from a binary image.
%
%   Input arguments:
%   - filename: ;bmp image that contain the bone structure.
%
%   Output arguments:
%   - boundaryEndost: Vector of the coordinates of the extracted boundary 
%   points. X-axis = (1, :) Y-axis = (2, :)
%   - boundaryPores: Cell array containing a double vector that contain the
%   X and Y coordinates of the pores.

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

    % PORES COMPUTATION
    boundaryPores = {};
    for i = 2:numel(boundaries)
        if numel(boundaries{i}) > 30
            boundaryPores{end+1} = flip(boundaries{i}, 2)';
        end
    end
end