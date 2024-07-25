function[Map, grid] = MakeGeometryExVivo(grid, probe, bone, rotate, withoutPores, printPlot, simu_dir)
% MakeGeometryExVivo generates a map from CT bone images.
%
% Syntax:
%   Map = MakeGeometryExVivo(grid, probe, bone, rotate, withoutPores, printPlot, simu_dir)
%
% Description:
%   MakeGeometryIExVivo generates a map from CT bone images.
%
% Input Arguments:
% Generated with GenerateAllParameters(). 
%   - grid: Structure containing information about the grid dimensions.
%   - probe: Structure containing information about the probe properties.
%   - bone : Structure containing the bone informations.
%   - rotate : Boleean to indicate is the bone needs to be rotated.
%   - printPlot: Logical value indicating whether to plot the excitation signal.
%   - simu_dir: Directory where the simulations files are saved, 
%           if not indicated, the signal is only computed on matlab, 
%           not registered in the simulation directory, the height profile is also saved in a .mat file.
%
% Output Arguments:
%   - Map: Binary map representing the bone/soft tissue interface.
%
% See also : imbinarize, SimSonic2DWriteMap2D
    
    to_px = @(mm) round(mm/grid.step); 
    % MAP GENERATION
    % Bone index :1
    % Soft tissu index : 0
    dirname = '~/Documents/BoneRugosity/2_ParametersDefinition/BoneImage';
    file = sprintf('SAMPLE_%s_SLICE_%04d.bmp', bone.id, bone.image); 
    filename = fullfile(dirname, file); 

    % Image binarization
    bone_bmp = imread(filename); 

    threshold = graythresh(bone_bmp); % Find an automatic threshold
    binaryImage = imbinarize(bone_bmp, threshold);
    
    % Rotate the image in order to obtain the endost on the bottom of the
    % map
    if rotate
        binaryImage = imrotate(binaryImage, -15,'bilinear','loose');
        binaryImage = binaryImage(1:find(...
            sum(binaryImage(: , 1: size(binaryImage, 2)/2),2) > 10, 1, 'last'), :);
    end 

    % The bone need to be in the middle of the image and with at least 2 mm
    % between the probe and the bone
    min_periost = find(sum(binaryImage, 2) > 0, 1);
    if min_periost < to_px(probe.depth+2)
        numLigneAdded = to_px(probe.depth+2) - min_periost;
    else 
        numLigneAdded = 0;
    end
    
    Nz = size(binaryImage, 1) + numLigneAdded; % Number of point in the direction 1 (depth - Z)
    Nx = to_px(grid.width);    % Number of point in the direction 2 (width - X)

    grid.depth = Nz * grid.step;

    Map = zeros(Nz, Nx, 'uint8');
    width = (Nx - size(binaryImage, 2))/2;

    if withoutPores
        binaryImage = imfill(binaryImage,'holes');
    end

    Map(numLigneAdded + 1 : end, width : end-width-1) = uint8(binaryImage);
    
    if nargin == 7
        fprintf(['\n--- Map and profile saved in ', simu_dir(46:end-1)]);
        SimSonic2DWriteMap2D(Map, fullfile(simu_dir, 'Geometry.map2D'));
    end
    
    % Plot map      
    if printPlot
        X = 0:grid.step:grid.width-grid.step; X = X -mean(X);
        Z = flipud(0:grid.step:grid.depth-grid.step);
        figure, imagesc(X, Z, Map)
        xticks(X(1:100:Nx));
        yticks(Z(1:100:Nz));
        axis equal
        xlabel('Width (mm)');
        ylabel('Depth (mm)');
        title('Simulation map');
    end
end


