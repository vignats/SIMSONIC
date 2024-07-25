function [surface, surfaceExtended] = ExpandParabola(boundaryEndost, boundaryWidth, option)
% EXPANDPARABOLA Expands the parabola fitted to the boundary points of the endosteum.
%
%   [surface, surfaceExtended] = ExpandParabola(boundaryEndost, boundaryWidth, option)
%   fits a parabolic curve to the boundary points of the endosteum obtained from the
%   binary image and then expands this parabola based on the specified boundary width.
%   It returns the coordinates of both the original and the expanded parabolas.
%
%   Inputs:
%   - boundaryEndost: Coordinates of the endosteum boundary.
%   - boundaryWidth: Width to expand the fitted parabola.
%   - option: Option to specify how the extended parabola should be computed.
%             'parameters' when the function is used to compute parameters, 
%             'specularity' when the function is used to determine 
%               the endosteal ROI for the specularity computation.
%
%   Outputs:
%   - surface: Matrix containing the coordinates of the fitted parabola.
%   - surfaceExtended: Matrix containing the coordinates of the expanded parabola.
%
%   Example:
%       [surface, surfaceExtended] = ExpandParabola(boundaryEndost, boundaryWidth, 'parameters')
%
%   See also: FitParabola 
    surface = FitParabola(boundaryEndost);
    
    % Get coefficient of the parabola
    % In order to compute he model, we reshape it. 
    y = max(surface(:,2)) - surface(:,2);
    x = surface(:,1) - unique(surface(y == max(y), 1));
    coeffParabola = polyfit(x, y, 2); 
    
    % Get maximum and root of the parabola
    r = max(roots(coeffParabola));
    
    % Multiply the coefficients to obtain the extanded parabola 
    A = (coeffParabola(1)*r^2 + boundaryWidth)/(r - boundaryWidth)^2;
    B = coeffParabola(2)*r/(r - boundaryWidth)^2;
    C = coeffParabola(3) - boundaryWidth;
    if coeffParabola(1) < 0
        A = coeffParabola(1);
        B = coeffParabola(2);
        C = coeffParabola(3) - boundaryWidth;
    end

    switch option
        case 'parameters'
            xExtended = (-length(x): 1: length(x))';
        case 'specularity'
            xExtended = x';
    end
    yExtended = A.*xExtended.^2 + B.*xExtended + C;

    % Reshape as previously 
    surfaceExtended(:,1) = xExtended + unique(surface(y == max(y), 1));
    surfaceExtended(:,2) = max(surface(:,2)) - yExtended - 2*boundaryWidth;

    surfaceExtended = surfaceExtended';
end