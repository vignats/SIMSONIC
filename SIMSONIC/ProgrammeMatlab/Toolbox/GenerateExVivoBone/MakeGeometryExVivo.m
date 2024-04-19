function[Map, heights] = MakeGeometryExVivo(grid, interface, filter, print, simu_dir)
% MakeGeometryWaviness generates a map representing the bone/soft tissue interface.
%
% Syntax:
%   Map = MakeGeometryWaviness(grid, probe, medium, interface, waviness, print, simu_dir)
%
% Description:
%   MakeGeometryInterface generates a map representing the bone/soft tissue interface based on extraction of the endost of the bone define in waviness.
%   The interface is obtained with GetRoughness, regarding the cut-off frequency.
%
% Input Arguments:
% Generated with GenerateAllParameters(). 
%   - grid: Structure containing information about the grid dimensions.
%   - interface: Structure containing information about the interface properties.
%   - waviness : Structure containing the bone information and cut-off
%   frequency. fc = 0.06 for waviness + roughness and 1.25 for roughness.
%   - print: Logical value indicating whether to plot the excitation signal.
%   - simu_dir: Directory where the simulations files are saved, 
%           if not indicated, the signal is only computed on matlab, 
%           not registered in the simulation directory, the height profile is also saved in a .mat file.
%
% Output Arguments:
%   - Map: Binary map representing the bone/soft tissue interface.
%
% See also : GetRoughness, SimSonic2DWriteMap2D

    to_px = @(mm) round(mm/grid.step); 

    % INTERFACE GENERATION
    [profile, heightInitial, xProfile] = GenerateInterface(filter);

    % Extend vector to recover the initial perimeters 
    heightExtended = ExtendInterface(profile, heightInitial, grid);

    % Make symmetry in order to obtain the desired width
    heights = MakeSymmetry(heightExtended, grid);

    % MAP GENERATION
    % Bone index :1
    % Soft tissu index : 0
    Nz = to_px(grid.depth);    % Number of point in the direction 1 (depth - Z)
    Nx = to_px(grid.width);    % Number of point in the direction 2 (width - X)
    Map=zeros(Nz, Nx, 'uint8');

    for i = 1:Nx
        Map(1: to_px(heights(i) + interface.depth), i) = 1;
    end 
    
    if nargin == 5
        fprintf(['\n--- Map and profile saved in ', simu_dir(46:end-1)]);
        SimSonic2DWriteMap2D(Map, [simu_dir 'Geometry.map2D']);
        save(fullfile(simu_dir, 'interface.mat'), "heightInitial", "profile", "heights");
    end
    
    % Plot map      
    if print
        X = 0:grid.step:grid.width-grid.step; X = X -mean(X);
        Z = flipud(0:grid.step:grid.depth-grid.step);
        figure, imagesc(X, Z, Map)
        xticks(X(1:100:Nx));
        yticks(Z(1:100:Nz));
        axis equal
        xlabel('Width (mm)');
        ylabel('Depth (mm)');
        title('Simulation map');

        figure
        plot(xProfile, heightInitial)
        title('Roughness of the interface, fc = ', filter.fc);
        xlabel('Width (mm)');
        ylabel('Depth (mm)');
        axis equal
        
        format = 'The interface have a rms of %.3f mm and a correlation length of %.3f mm \n';
        fprintf(format, rms(heights), ComputeCorr(heights, grid.step));
    end
end

function [profile, heightInitial, xProfile] = GenerateInterface(filter)
% This function allows to generate a planar interface with rugosity, based
% on an ex-vivo bone. 
    dirname = ['/calculSSD//Dossier partagé image os exvivo/', filter.bone, '/'];
    if filter.segmented
        dirname = [dirname, 'SEGMENTED_OTSU3D_K4/'];
    end
    file = ['SAMPLE_', filter.bone, '_SLICE_', sprintf('%04d', filter.image), '.bmp']; 
    filename = [dirname, file];
    
    [profile, heightInitial, xProfile] = FilterProfil(filename, filter.fc, filter.segmented);
end

function [heightExtended] = ExtendInterface(profile, heightInitial, grid)
% This function applies a coefficient to the filtered profile in order to 
% recover the initial distances between the points of the initial curved 
% profile. 
    perimetreProfile = sum(sqrt(grid.step.^2 + diff(profile).^2));
    xExtended = (1:perimetreProfile/grid.step)*grid.step;
    xInitial = (1:length(heightInitial))*perimetreProfile/length(heightInitial);
    heightExtended = interp1(xInitial,heightInitial,xExtended);
    heightExtended(1) = heightInitial(1);
end

function [heights] = MakeSymmetry(heightExtended, grid)
% This function allow to add the symmetrie of the profile so that the
% width of the profile is equal to the grid width.
    to_px = @(mm) round(mm/grid.step); 

    % Compute symmetrie and repeat it to be above the gird width
    % to be accurate we should use a step of filter.fs = 9e-3mm, we are
    % here introducing an error of 1um at each pixel.
    nbRepetiton = ceil(grid.width/(2*length(heightExtended)*grid.step));
    heightSymmetry = [fliplr(heightExtended), heightExtended];
    heightWiden = repmat(heightSymmetry, 1, nbRepetiton);

    heights = heightWiden(1: to_px(grid.width));
end

