function [porosity, poreSize, imageDisplay] = ComputePorosity(binaryImage, boundaryEndost, boundaryPores, boundaryWidth, plotBoundary)
% This function allows to compute a percentage of porosity which is
% represented by the volume of bone over the total volume in a specified
% width around the boundary.
% 
% Input :   binaryImage - Binarized image of an ex-vivo X-Rayed image.
%           boundaryEndost - Boundary of the endost, usually extracted with
%           the function ExtractBoundary.
%           boundaryWidth - Width of the boundary in which the porosity is
%           computed, in mm.
%           plotBoundary - boolean that indicate weither to plot the
%           extended boundary.
% 
% Output :  porosity - Percentage of bone in the total volume = BV/TV *100
%
% See also : ExpandParabola
    
    % Define the width corresponding to a number of wavelength
    pixelSize = 9e-3; 
    boundaryWidth = boundaryWidth/pixelSize;
    
    [surface, surfaceExtended] = ExpandParabola(boundaryEndost, boundaryWidth, 'parameters');
    
    boneVolume = 0;
    totalVolume = 0;
    imageDisplay = double(binaryImage);
    
    % Compute Pore Size
    poresIn = {};
    surfaceIn = surfaceExtended(1, surfaceExtended(2,:) < boundaryEndost(2,1));
    
    for j = 1 : length(boundaryPores)
        xMean = mean(boundaryPores{j}(1,:));
        yMean = mean(boundaryPores{j}(2,:));
        
        if xMean > surfaceIn(1) && xMean < surfaceIn(end) 
            % Trouver les indices correspondants pour xMean
            indexSurface = find(surfaceExtended(1,:) >= xMean, 1);
            indexBoundary = find(boundaryEndost(1,:) >= xMean, 1);
            
            % S'assurer que les indices sont valides
            if ~isempty(indexSurface) && ~isempty(indexBoundary)
                ySurface = surfaceExtended(2, indexSurface);
                yBoundary = boundaryEndost(2, indexBoundary);
                
                if yMean > ySurface && yMean < yBoundary
                    poresIn{end+1} = boundaryPores{j};
                end
            end
        end
    end
    
    poreLength = zeros(length(poresIn),1);
    for i = 1:length(poresIn)
        poreLength(i) = length(poresIn{i});
    end
    poreSize = mean(poreLength)*pixelSize/pi;
    
    % Compute porosity
    for i = surfaceExtended(1, surfaceExtended(2,:) < boundaryEndost(2,1))
        if i < boundaryEndost(1,1) || i > boundaryEndost(1,end)
            limSup = max(boundaryEndost(2,:));
        else 
            limSup = max(boundaryEndost(2,boundaryEndost(1,:) == i));
        end
        
        boneVolume = boneVolume + sum(binaryImage(ceil(surfaceExtended(2,surfaceExtended(1,:) == i)) : limSup, i));
        totalVolume = totalVolume + numel(binaryImage(ceil(surfaceExtended(2,surfaceExtended(1,:) == i)) : limSup, i));

        % Change the value of the pixel in the boundary 
        imageDisplay(ceil(surfaceExtended(2,surfaceExtended(1,:) == i)) : limSup, i) =...
            imageDisplay(ceil(surfaceExtended(2,surfaceExtended(1,:) == i)) : limSup, i).*0.5;
    end

    porosity = 100 * (1 - boneVolume/totalVolume);

    if plotBoundary
        figure
        imshow(binaryImage);
        hold on 
        plot(boundaryEndost(1,:), boundaryEndost(2,:), 'r', 'LineWidth', 2)
        hold on         
        plot(surfaceExtended(1,:), surfaceExtended(2,:), 'b', 'LineWidth', 2)
        hold on  
        plot(surface(:,1), surface(:,2), 'g', 'LineWidth', 1)
        for k = 1:length(poresIn)
            hold on
            plot(poresIn{k}(1,:), poresIn{k}(2,:), 'LineWidth', 2);
        end
        
        xlabel('Width', 'Interpreter', 'latex', 'FontSize', 22);
        ylabel('Depth', 'Interpreter', 'latex', 'FontSize', 22);
        grid on
        title('Boundary of the endostium', 'Interpreter', 'latex', 'FontSize', 24); 
        legend('Boundary of the endost', 'Expended boundary', 'Fitted Parabola', 'Interpreter', 'latex', 'FontSize', 18); 
        hold off;
        ax = gca; 
        ax.FontSize = 16;
        
    end
end
%% Script to plot the expanded boundary. 
% y = max(surface(:,2)) - surface(:,2);
% x = surface(:,1) - unique(surface(y == max(y), 1));
% angle = atan2(surface(end,2) - surface(1,2), surface(end,1) - surface(1,1));
% 
% % Matrice de rotation
% R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
% 
% % Appliquer la rotation à tous les points de la parabole
% rotated_points = R * [x'; y'];
% 
% % Séparer les points x et y
% x_rotated = rotated_points(1, :);
% y_rotated = rotated_points(2, :);
% 
% x_rotated = x_rotated + unique(surface(y == max(y), 1));
% y_rotated = max(surface(:,2)) - y_rotated;
% 
% x = x + unique(surface(y == max(y), 1));
% y = max(surface(:,2)) - y;
% 
% % Ajuster pour que le début et la fin soient à zéro
% figure;
% plot(x, y, 'b-', 'DisplayName', 'Original Parabola'); % Parabole originale
% hold on;
% plot(x_rotated, y_rotated, 'r-', 'DisplayName', 'Rotated Parabola'); % Parabole pivotée
% legend show;
% xlabel('x');
% ylabel('y');
% title('Rotation de la Parabole');
% grid on;
% 
% figure
% imshow(binaryImage);
% hold on 
% plot(boundaryEndost(1,:), boundaryEndost(2,:), 'r', 'LineWidth', 2)
% hold on         
% plot(surface(:, 1), surface(:,2), 'b', 'LineWidth', 2)
% hold on         
% plot(surface(:, 1), y_rotated)

%%

% for i = surfaceExtended(1, surfaceExtended(2,:) < boundaryEndost(2,1))
%     if i < boundaryEndost(1,1) || i > boundaryEndost(1,end)
%         limSup = max(boundaryEndost(2,:));
%     else 
%         limSup = max(boundaryEndost(2,boundaryEndost(1,:) == i));
%     end
% 
%     poresIn = {};
%     for j = 1 : length(boundaryPores)
%         xMean = mean(boundaryPores{j}(1,:));
%         yMean = mean(boundaryPores{j}(2,:));
%         if xMean > a(1) && xMean > a(end) 
%             if yMean > surfaceExtended(find(surfaceExtended(1,:) >= xMean, 1), 2) &&...
%                     yMean < boundaryEndost(find(boundaryEndost(1,:) >= xMean, 1), 2)
%                 poresIn{end+1} = boundaryPores{j};
%             end
%         end
%     end
% 
% end
% figure
% imshow(binaryImage);
% hold on 
% plot(boundaryEndost(1,:), boundaryEndost(2,:), 'r', 'LineWidth', 2)
% hold on         
% plot(surfaceExtended(1,:), surfaceExtended(2,:), 'b', 'LineWidth', 2)
% hold on         
% plot(surface(:,1), surface(:,2), 'g', 'LineWidth', 2)
% hold on
% for k = 1:length(boundaryPores)
%     plot(boundaryPores{k}(1,:), boundaryPores{k}(2,:));
%     hold on
% end
% xlabel('Width');
% ylabel('Depth'); 