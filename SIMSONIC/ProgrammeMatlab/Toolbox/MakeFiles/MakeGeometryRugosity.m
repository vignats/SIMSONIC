function[Map] = MakeGeometryRugosity(grid, probe, medium, interface, print_plot, print_hist, simu_dir)
	to_px = @(mm) round(mm/grid.step); 
    Nz = to_px(grid.depth);    % Number of point in the direction 1 (depth - Z)
    Nx = to_px(grid.width);    % Number of point in the direction 2 (width - X)
    
    % MAP WITHOUT POROSITY
    % Bone index :1
    % Soft tissu index : 0
    Map=zeros(Nz, Nx, 'uint8');
    Map(1 : to_px(interface.endost), :) = 1;
    
    % PORE GENERATION
    % The rugosity corresponds to the volume of pore in a wavelength width after 
    % bone/soft tissu interface 
    lambda = (medium.cp(2) / (probe.fc))*1e3;      % Wavelength in the bone (mm)
    bone_area = Nx*to_px(lambda);                   % In pixel 
    total_pore_area = bone_area * interface.rugosity/100; 
    pore_area = 0;
    
    % Distribution 
    position_dist = [];
    size_dist = [];
    [x, z] = meshgrid(1:Nx, 1:Nz);
    tic
    h = waitbar(0, 'Geometry map processing...');

    while pore_area <= total_pore_area
        [pore_x, pore_z, pore_size] = MakePore(grid, interface);
        
        position_dist(end+1) = pore_z;
        size_dist(end+1) = pore_size;
        
        % Modification of the Map to add pore
        pore_x_px = to_px(pore_x);
        pore_z_px = to_px(pore_z);
        pore_size_px = to_px(pore_size/2);
        
        masque = (x - pore_x_px).^2 + (z - pore_z_px).^2 >= pore_size_px.^2;            % Création of a mask to erode
        Map = logical(Map).* masque;   % max and min to make sure we are not out of the Map index
   
        pore_area = sum(Map(to_px(interface.endost) - to_px(lambda) :to_px(interface.endost), :) == 0, 'all');
        
        progress = pore_area / total_pore_area;
        waitbar(progress, h, sprintf('Progress: %.f%%', progress*100));
    end 
    close(h);
    toc
    
    if nargin == 7
        SimSonic2DWriteMap2D(uint8(Map),[simu_dir 'Geometry.map2D']);
    end
    
    % Plot map      
    if print_plot
        X = 0:grid.step:grid.width-grid.step; X = X - mean(X);
        Z = flipud(0:grid.step:grid.depth-grid.step);
        figure, imagesc(X, Z, Map)
        xticks(X(1:100:Nx));
        yticks(Z(1:100:Nz));
        axis equal
        xlabel('Width (mm)');
        ylabel('Depth (mm)');
        title('Simulation map');
        
        % Rugosity calculation 
        [Rq, Corr] = ComputeParam(Map, grid);
        rugosity = ComputeRugosity(Map, grid, probe, medium, interface);
        format = 'The interface have a rms of %.3f mm and a correlation length of %.3f mm \nThe rugosity in one wavelength is %.3f';
        fprintf(format,Rq, Corr, rugosity);
    end
    
    % Plot histogram of simulated diameters
    if print_hist
        figure;
        histogram(size_dist);
        xlabel('Diameter (mm)');
        ylabel('Number of pores');
        title('Simulated Diameter Distribution');
        
        figure;
        histogram(position_dist);
        xlabel('Position');
        ylabel('Number of pores');
        title('Simulated Position Distribution');
    end 
end

function [pore_x, pore_z, pore_size] = MakePore(grid, interface)
    % This function return the coordinate of a pore and its diameter. 
    
    % PORE SIZE 
    % The pore size follow a gamma distribution, the alpha and beta
    % parameters depends of the patients. 

    if strcmp(interface.patient, 'aged')
        alpha = 25; beta = 2;                       
    elseif strcmp(interface.patient, 'osteoporotic')
        alpha = 10; beta = 7;    
    elseif strcmp(interface.patient, 'young')
        alpha = 35;  beta = 1; 
    else 
        disp('Unknown patient type');
    end
                         
    pore_size = gamrnd(alpha, beta)*1e-3;     % Alpha and beta values corresponds to values in um        

    % PORE POSITION
    % Position of the pore along the X-axis 
    pore_x = rand*grid.width;

    % Position of the pore is gaussian distributed along the Z-axis,
    % centered at the medium interface
    mu = interface.endost;                           % Center of the distribution (mm)
    std = 0.5*1.1750;                               % FWHM/2 = 0.5 mm for the pore distribution and FWHM = 2.35*std 

    pore_z = mu - abs(normrnd(0, std));
end 

function [rugosity] = ComputeRugosity(Map, grid, probe, medium, interface)
   % The rugosity corresponds to the volume of pore in the 200um width after 
   % bone/soft tissu interface 
   to_px = @(mm) round(mm/grid.step);
   lambda = (medium.cp(2) / (probe.fc))*1e3;       % Wavelength in the bone (mm)
   lambda = to_px(lambda);
   
   pore_area = sum(Map(to_px(interface.endost) - lambda :to_px(interface.endost), :) == 0, 'all');
   tot_area = numel(Map(to_px(interface.endost) - lambda :to_px(interface.endost), :));
   
   rugosity = pore_area / tot_area *100;
end

function [Rq, Corr] = ComputeParam(Map, grid)
% This function compute the root mean square and correlation length that
% describe the interface. 
    % Extract the bone interface
    profile = zeros(1, size(Map, 2));
    for i = 1:size(Map, 2)
        profile(1, i) = find(Map(:, i) == 1, 1, 'last'); % Find last bone element that correspond to the bone ST interface
    end
    
    Rq = rms(profile)*grid.step;   % RMS of the surface (mm)   
    
    Corr = ComputeCorr(profile, grid.step);
end

