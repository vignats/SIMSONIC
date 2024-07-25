function[Map, heights] = MakeGeometry(grid, probe, medium, interface, verify, simu_dir)
% MakeGeometry generates a map representing the bone/soft tissue interface.
%
% Syntax:
%   Map = MakeGeometry(grid, probe, medium, interface, simu_dir, print)
%
% Description:
%   MakeGeometry generates a map representing the bone/soft tissue interface based on the provided parameters.
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
    iteration = 0; % Set a limit to the number of iteration. 
    
    while regenerate
    
        % MAP WITHOUT POROSITY
        % Bone index :1
        % Soft tissu index : 0
        Map=zeros(Nz, Nx, 'uint8');
        
        heights = GenerateRoughSurface(Nx, grid.width, interface.rms, interface.corr);
        heights(isnan(heights)) = 0;
        
        periInd = 1;
        if interface.periost ~= 0
            periInd = to_px(interface.periost);
        end 
        for i = 1:Nx
            Map(periInd: to_px(heights(i) + interface.endost), i) = 1;
        end 
        
        % Verify that the plot verify the rms and correlation length asked
        [Rq, Corr, ~] = ComputeInterfaceParameters(Map, grid, probe, medium);
        iteration = iteration +1;
        if (0.9*interface.rms < Rq) && (Rq < 1.1*interface.rms) ...
                && (0.9*interface.corr < Corr) && (Corr < 1.1*interface.corr) ...
                || iteration >100 || interface.rugosity == 0
            regenerate = false;
            if iteration > 100
                frpintf('Not able to create map with input parameters')
            end
        end
        
        % Generate pore inside the bone
        if interface.rugosity > 0
            % Compute the number of pore needed in a width of one
            % wavelength above the bone/ST interface. 
            lambda = medium.cp(2)*1e3 / (probe.fc);     % Wavelength in the bone (mm)
            
            maxRowBone = find(Map(:, find(sum(Map)==max(sum(Map)), 1)) == 1, 1, 'last'); 
            totalVolume = sum(Map(to_px(interface.endost - lambda) : maxRowBone, :), 'all');

            % nbPore = round((totalVolume*interface.rugosity/100)/(pi*to_px(interface.diameter)^2/4)); % sum(Map(to_px(interface.endost) - lambda : to_px(interface.endost)))
            % xPores = randi(Nx, nbPore, 1);   
            % zPores = randi([to_px(interface.endost - lambda), maxRowBone], nbPore, 1); 
            % for i = 1:nbPore
            
            % Generate pore until we reached the desired pore volume
            [x, z] = meshgrid(1:Nx, 1:Nz);
            while sum(Map(to_px(interface.endost - lambda) : maxRowBone, :), 'all') >= totalVolume*(1-interface.rugosity/100)
                % Generate a X and Z-coordinate for the pore between the
                % original interface minus a wavelength and the bone endost
                % generated previously.
                xPores = randi(Nx, 1, 1);   
                zPores = randi([to_px(interface.endost - lambda), maxRowBone], 1, 1);

                % Once a unique paire of coordinate has been generated,
                % transform the bone in soft tissu medium (0) in a specified radius around the pore.
                masque = (x - xPores).^2 + (z - zPores).^2 >= to_px(interface.diameter/2).^2;            % Création of a mask to erode
                Map = logical(Map).* masque;
                % Map(zPores - 1 : zPores + 1 , xPores - 1 :  xPores + 1) = 0;   
            end
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
        title('Simulation map', ...
            sprintf('rms = %.2f mm, corr = %.2f mm, rugosity = %.0f %%, pore diameter = %.0f um', ...
            interface.rms, interface.corr, interface.rugosity, interface.diameter*1e3));
        
        % Rugosity calculation 
        [Rq, Corr, rugosity] = ComputeInterfaceParameters(Map, grid, probe, medium);
        format = ['\nThe interface have a rms of %.3f mm and a correlation length of %.3f mm' ...
            '\nThe rugosity in one wavelength is %.1f%%'];
        fprintf(format, Rq, Corr, rugosity);
    end 

    % Compute Geometry.map2D file
    if nargin == 6
        fprintf(['\n--- Geometry Map saved in ', simu_dir(46:end-1)]);
        SimSonic2DWriteMap2D(uint8(Map), fullfile(simu_dir, 'Geometry.map2D'));
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