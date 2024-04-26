function [porosity, totalVolume] = ComputePorosity(binaryImage, boundaryEndost, nbWavelength, plotBoundary)
% This function allows to compute a percentage of porosity which is
% represented by the volume of bone over the total volume in a specified
% width around the boundary.
% 
% Input :   binaryImage - Binarized image of an ex-vivo X-Rayed image.
%           boundaryEndost - Boundary of the endost, usually extracted with
%           the function ExtractBoundary.
%           nbWavelength - Width of the boundary in which the porosity is
%           computed, in number of wavelength inside the bone.
%           plotBoundary - boolean that indicate weither to plot the
%           extended boundary.
% 
% Output :  porosity - Percentage of bone in the total volume = BV/TV *100
%
% See also : ExpandParabola
    
    % Define the width corresponding to a number of wavelength
    pixelSize = 9e-3; lambda = 1.4; % Wavelength of the ultrasound in the bone (mm) 
    boundaryWidth = nbWavelength*lambda/pixelSize;

    surfaceExtended = ExpandParabola(boundaryEndost, boundaryWidth);

    boneVolume = 0;
    totalVolume = 0;

    for i = surfaceExtended(1, surfaceExtended(2,:) < boundaryEndost(2,1))
        if i < boundaryEndost(1,1) || i > boundaryEndost(1,end)
            limSup = max(boundaryEndost(2,:));
        else 
            limSup = boundaryEndost(2,boundaryEndost(1,:) == i);
        end
        
        boneVolume = boneVolume + sum(binaryImage(ceil(surfaceExtended(2,surfaceExtended(1,:) == i)) : limSup, i));
        totalVolume = totalVolume + numel(binaryImage(ceil(surfaceExtended(2,surfaceExtended(1,:) == i)) : limSup, i));
    end

    porosity = 100 * (1 - boneVolume/totalVolume);

    if plotBoundary
        figure
        imshow(binaryImage);
        hold on 
        plot(boundaryEndost(1,:), boundaryEndost(2,:), 'r', 'LineWidth', 2)
        hold on         
        plot(surfaceExtended(1,:), surfaceExtended(2,:), 'b', 'LineWidth', 2)
        xlabel('Width');
        ylabel('Depth');
        grid on
        title('Boundary of the endostium'); 
        legend('Boundary of the endost', 'Limite of the boundary'); hold off;

        fprintf('The bone as a porosity of %.1f %% in a width of %.1f wavelength before the endost thus %.2f mm', ...
            porosity, nbWavelength, nbWavelength*lambda);
    end
end

function [surfaceExtended] = ExpandParabola(boundaryEndost, boundaryWidth)
% This function allows to expand the parabola that fit to the boundary
% endost. 
    surface = FitParabola(boundaryEndost);

    % Get coefficient of the parabola
    % In order to compute he model, we reshape it. 
    y = max(surface(:,2)) - surface(:,2);
    x = surface(:,1) - surface(y == max(y), 1);
    coeffParabola = polyfit(x, y, 2); 
    
    % Get maximum and root of the parabola
    r = max(roots(coeffParabola));
    
    % Multiply the coefficients to obtain the extanded parabola 
    A = (coeffParabola(1)*r^2 + boundaryWidth)/(r - boundaryWidth)^2;
    B = coeffParabola(2)*r/(r - boundaryWidth)^2;
    C = coeffParabola(3) - boundaryWidth;
    
    xExtended = (-length(x): 1: length(x))';
    yExtended = A.*xExtended.^2 + B.*xExtended + C;

    % Reshcape as previously 
    surfaceExtended(:,1) = xExtended + surface(y == max(y), 1);
    surfaceExtended(:,2) = max(surface(:,2)) - yExtended - 2*boundaryWidth;

    surfaceExtended = surfaceExtended';
end
%% Script to plot the expanded boundary. 
