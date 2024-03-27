function[Map] = MakeGeometryWaviness(grid, probe, medium, interface, filter, print, simu_dir)
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
%   - probe: Structure containing information about the ultrasound probe.
%   - medium: Structure containing information about the medium properties.
%   - interface: Structure containing information about the interface properties.
%   - waviness : Structure containing the bone information and cut-off
%   frequency. fc = 0.06 for waviness + roughness and 1.25 for roughness.
%   - print: Logical value indicating whether to plot the excitation signal.
%   - simu_dir: Directory where the simulations files are saved, 
%           if not indicated, the signal is only computed on matlab, 
%           not registered in the simulation directory.
%
% Output Arguments:
%   - Map: Binary map representing the bone/soft tissue interface.

    to_px = @(mm) round(mm/grid.step); 
    Nz = to_px(grid.depth);    % Number of point in the direction 1 (depth - Z)
    Nx = to_px(grid.width);    % Number of point in the direction 2 (width - X)

    % MAP WITHOUT POROSITY
    % Bone index :1
    % Soft tissu index : 0
    Map=zeros(Nz, Nx, 'uint8');

    % INTERFACE GENERATION
    % heigths = GenerateRoughSurface(Nx, grid.width, interface.rms, interface.corr);
    heigths = GenerateInterface(filter);

    for i = 1:Nx
        Map(1: to_px(heigths(i) + interface.depth), i) = 1;
    end 
    
    if nargin == 6
        SimSonic2DWriteMap2D(Map, [simu_dir 'Geometry.map2D']);
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
        
        % Rugosity calculation 
        format = 'The interface have a rugosity of %.1f %%\n';
        fprintf(format, ComputeRugosity(Map, grid, probe, medium, interface));
    end
end

function [heigths] = GenerateInterface(filter)
    dirname = ['~/Documents/BoneRugosity/RMS/', filter.bone, '/'];
    file = ['SAMPLE_', filter.bone, '_SLICE_', num2str(filter.image), '.bmp']; 
    filename = [dirname, file];
    
    % Load image
    bone_bmp = imread(filename); 

    heigths = GetRoughness(bone_bmp, filter.fc);
end