function[Map] = MakeGeometryInterface(grid, probe, medium, interface, print, simu_dir)
% MakeGeometryInterface generates a map representing the bone/soft tissue interface.
%
% Syntax:
%   Map = MakeGeometryInterface(grid, probe, medium, interface, simu_dir, print)
%
% Description:
%   MakeGeometryInterface generates a map representing the bone/soft tissue interface based on the provided parameters.
%   The interface is obtained with rsdene1D, regarding the RMS and correlation length.
%
% Input Arguments:
% Generated with GenerateAllParameters(). 
%   - grid: Structure containing information about the grid dimensions.
%   - probe: Structure containing information about the ultrasound probe.
%   - medium: Structure containing information about the medium properties.
%   - interface: Structure containing information about the interface properties.
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
    heigths = GenerateRoughSurface(Nx, grid.width, interface.rms, interface.corr);

    for i = 1:Nx
        Map(1: to_px(heigths(i) + interface.depth), i) = 1;
    end 
    
    if nargin == 6
        SimSonic2DWriteMap2D(Map, [simu_dir 'Geometry.map2D']);
    end
    
    % Plot map      
    if print
        X = 0:grid.step:grid.width-grid.step;
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

function [heights] = GenerateRoughSurface(N, length, rms, corr)
% From rsgene1D.m function by David Bergström
% Last updated: 2010-07-26 ()
%
% [heights] = GenerateRoughSurface(N,grid, interface)
%
% generates a 1-dimensional random rough surface f(x) with N surface points.
% The surface has a Gaussian height distribution and an
% exponential autocovariance function, where rL is the length of the surface, 
% h is the RMS height and cl is the correlation length.
%
% Input:    N       - number of surface points
%           length  - length of the surface 
%           rms     - RMS of the height distribution 
%           corr    - correlation length
%
% Output:   heights   - surface heights

    format long;
    
    x = linspace(-length/2,length/2,N);
    
    Z = rms.*randn(1,N); % uncorrelated Gaussian random rough surface distribution
                       % with rms height h
                            
    % Gaussian filter
    F = exp(-abs(x)/(corr/2));
    
    % correlated surface generation including convolution (faltning) and inverse
    % Fourier transform and normalizing prefactors
    heights = sqrt(2)*sqrt(length/N/corr)*ifft(fft(Z).*fft(F));
end

function [rugosity] = ComputeRugosity(Map, grid, probe, medium, interface)
   % The rugosity corresponds to the volume of pore in the 200um width after 
   % bone/soft tissu interface 
   to_px = @(mm) round(mm/grid.step);
   lambda = (medium.cp(2) / (probe.fc*1e3));       % Wavelength in the bone (mm)
   lambda = to_px(lambda);
   
   pore_area = sum(Map(to_px(interface.depth) - lambda :to_px(interface.depth), :) == 0, 'all');
   tot_area = numel(Map(to_px(interface.depth) - lambda :to_px(interface.depth), :));
   
   rugosity = pore_area / tot_area *100;
end