function[boundaryEndost, boundaryPores] = ExtractBoundary(binaryImage)
%   This function extracts the boundary of the endost from a binary image.
%
%   Input arguments:
%   - binaryImage: binary image that contain the bone structure.
%
%   Output arguments:
%   - boundaryEndost: Vector of the coordinates of the extracted boundary 
%   points. X-axis = (1, :) Y-axis = (2, :)
%   - boundaryPores: Cell array containing a double vector that contain the
%   X and Y coordinates of the pores.

    % Select boundary of the endost and periost  
    boundaries = bwboundaries(binaryImage, 8);
    if length(boundaries{1}) < 1000
        boundary = boundaries{1};
        for idx = 1 : length(boundaries)
            if length(boundaries{idx}) > length(boundary)
                boundary = boundaries{idx};
            end
        end
        boundary =  flip(boundary, 2);
    else
        boundary = flip(boundaries{1}, 2);      % First column X-axis & second column Z-axis
    end
    
    % BOUNDARIES COMPUTATION
    % Delimitation of the endost
    % Deletion of the periost by selecting the maximum value along each column.
    [~, idx] = unique(boundary(:,1), 'last');        % Get boundary index 
    X_bound = flip(boundary(min(idx) : max(idx), 1));      % Get X value for each column. 
    Z_bound = flip(boundary(min(idx) : max(idx), 2));      % Get corresponfing Z value boundary(idx, 2) and create a matrix of the endost coordinate

    % Deletion of the image boundary
    [nbOccurences, Zcoordinate] = groupcounts(Z_bound);                 % Count the number of occurences to detect the image boundary
    ImageBoundary = Zcoordinate(nbOccurences == max(nbOccurences));     % Define the Z coordinate of the image boundary
    
    Limite = (X_bound(Z_bound == ImageBoundary));                       % Find the X coordinates corresponding to the limite 
    LimEndoMax =  find(X_bound == min(Limite(Limite - mean(X_bound) > 0)), 1, 'last');   
    LimEndoMin =  find(X_bound == max(Limite(Limite - mean(X_bound) < 0)), 1);

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