function[Map, heights] = MakeGeometryInterface(grid, probe, medium, interface, verify, simu_dir)
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

    % INTERFACE GENERATION
    % While both calculated rms and correlation length aren't in a interval
    % of 10% the demanded value, we regenerate the Map.
    regenerate = true;
    intervalRms = [interface.rms - 0.1*interface.rms, interface.rms + 0.1*interface.rms];
    intervalCorr = [interface.corr - 0.1*interface.corr, interface.corr + 0.1*interface.corr];
    iteration = 0; % Set a limit to the number of iteration. 
    
    while regenerate
    
        % MAP WITHOUT POROSITY
        % Bone index :1
        % Soft tissu index : 0
        Map=zeros(Nz, Nx, 'uint8');

        heights = GenerateRoughSurface(Nx, grid.width, interface.rms, interface.corr);
        heights(isnan(heights)) = 0;
    
        for i = 1:Nx
            endostInd = 1;
            if interface.periost ~= 0
                endostInd = to_px(interface.periost);
            end 
            Map(endostInd: to_px(heights(i) + interface.endost), i) = 1;
        end 
        
        % Verify that the plot verify the rms and correlation length asked
        [Rq, Corr, ~] = ComputeInterfaceParameters(Map, grid, probe, medium);
        iteration = iteration +1;
        if (intervalRms(1) < Rq) && (Rq < intervalRms(2)) ...
                && (intervalCorr(1) < Corr) && (Corr < intervalCorr(2)) ...
                || iteration >100
            regenerate = false;
            if iteration > 100
                frpintf('Not able to create map with input parameters')
            end
        end
        % Plot map      
        if strcmp(verify, 'plot')
            X = 0:grid.step:grid.width-grid.step; X = X -mean(X);
            Z = flipud(0:grid.step:grid.depth-grid.step);
            figure, imagesc(X, Z, Map)
            xticks(X(1:100:Nx));
            yticks(Z(1:100:Nz));
            % axis equal
            xlabel('Width (mm)');
            ylabel('Depth (mm)');
            title('Simulation map', sprintf('rms = %.2f mm, corr = %.2f mm', interface.rms, interface.corr));
            
            % Rugosity calculation 
            [Rq, Corr, rugosity] = ComputeInterfaceParameters(Map, grid, probe, medium);
            format = '\nThe interface have a rms of %.3f mm and a correlation length of %.3f mm \nThe rugosity in one wavelength is %.1f%%';
            fprintf(format, Rq, Corr, rugosity);
        end
    end
        
    % Compute Geometry.map2D file
    if nargin == 6
        fprintf(['\n--- Geometry Map saved in ', simu_dir(46:end-1)]);
        SimSonic2DWriteMap2D(Map, [simu_dir 'Geometry.map2D']);
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
   
   pore_area = sum(Map(to_px(interface.endost) - lambda :to_px(interface.endost), :) == 0, 'all');
   tot_area = numel(Map(to_px(interface.endost) - lambda :to_px(interface.endost), :));
   
   rugosity = pore_area / tot_area *100;
end