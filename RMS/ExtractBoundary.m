function [boundary_endost, xy] = ExtractBoundary(binary_image, xy)
%   This function extracts the boundary of the endost from a binary image.
%   If the function is applied to the first image of the bone, the second
%   arguments need to be true, in this case : 
%       The function opens an interactive image to pre-select manually
%       the point that delimits the endost. 
%       Double clic on the left and right limite of the endost to selecte it.
%
%   Input arguments:
%   - BINARY_IMAGE: Binary image containing the bone structure.
%   - xy : matrix containing the coordinates of the points that define the 
%   limits of the endosteum. If the function is used on the first image, 
%   this argument should not be provided, allowing the function to compute 
%   these limits automatically.
%
%   Output arguments:
%   - BOUNDARY: Coordinates of the extracted boundary points.
%   - xy : matrix containing the coordinates of the points that define the 
%   limits of the endosteum, common to all the images of the same bone.
%
%   See also: getpts

    % Select boundary of the endost and periost    
    boundaries = bwboundaries(binary_image,8,'noholes');
    [~, max_index] = max(cellfun('size', boundaries, 1));
    boundary = flip(boundaries{max_index}, 2);
    
    % Delimitation of the endost using getpts, need to be done only once per
    % bone, for the first image
    if nargin == 1
        figure; 
        imshow(binary_image);
        hold on;
        plot(boundary(:, 1), boundary(:,2), 'o');

        [x1,y1] = getpts; 
        [x2,y2] = getpts; 
    
        xy = [x1 y1; x2 y2];
    end

    % The limitation of the endost need to be founded for each image of the
    % bone
    [lim, ~] = dsearchn(boundary,xy);   % Search the closest point to the manually selected limit of the endost
    boundary_endost = boundary(min(lim): max(lim), :);

    % Plot to verify for the first image
    if nargin == 1
        imshow(binary_image);
        hold on;
        plot(boundary_endost(:, 1), boundary_endost(:, 2), 'r', 'LineWidth', 2);
    end
end