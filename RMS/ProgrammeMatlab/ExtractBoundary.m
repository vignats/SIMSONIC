function [boundary_endost] = ExtractBoundary(binary_image)
%   This function extracts the boundary of the endost from a binary image.
%
%   Input arguments:
%   - BINARY_IMAGE: Binary image containing the bone structure.
%
%   Output arguments:
%   - BOUNDARY: Coordinates of the extracted boundary points.
%
%   See also: getpts

    % Select boundary of the endost and periost    
    boundaries = bwboundaries(binary_image, 8, 'noholes');
    [~, max_index] = max(cellfun('size', boundaries, 1));
    boundary = flip(boundaries{max_index}, 2);
    
    % BOUNDARIES COMPUTATION
    % Delimitation of the endost
    % Deletion of the periost.
    [X_bound, idx] = unique(boundary(:,1), 'last');
    boundary_endost = [X_bound, boundary(idx, 2)];

    boundary_endost = removerows(boundary_endost, 'ind', boundary_endost(:,2) > 660);
    boundary_endost = removerows(boundary_endost, 'ind', find(diff(boundary_endost(:,1)) > 100));
    boundary_endost = boundary_endost';
end