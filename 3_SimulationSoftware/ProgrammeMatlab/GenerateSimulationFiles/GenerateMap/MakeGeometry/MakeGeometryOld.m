function[] = MakeGeometryRealistic(grid, medium, pore, simu_dir, print, print_hist)
    to_px = @(mm) round(mm/grid.step); 
    Nz = to_px(grid.depth);    % Number of point in the direction 1 (depth - Z)
    Nx = to_px(grid.width);    % Number of point in the direction 2 (width - X)

    % Bone index :1
    % Soft tissu index : 0
    Map=zeros(Nz, Nx, 'uint8');
    Map(1 : to_px(medium.interface), :) = 1;

    % Pore generation
    [pore_x_px, pore_z_px, pore_size_px] = MakePorePorosity(grid, medium, pore, print_hist);
    for i = 1:length(pore_size_px)
        % Bone index to soft tissu index
        x = pore_x_px(i);
        z = pore_z_px(i);
        r = pore_size_px(i);  % The pore size corresponds to the diameter

        Map(z-r : z+r, max(1,x-r) : min(Nx, x+r)) = 0;  % max and min to make sure we are not out of the Map index
    end 

    SimSonic2DWriteMap2D(Map,[simu_dir 'Geometry.map2D']);

    % Plot map 
    if print == true
        X = 0:grid.step:grid.width-grid.step;
        Z = flipud(0:grid.step:grid.depth-grid.step);
        figure, imagesc(X, Z, Map)
        xticks(X(1:100:Nx));
        yticks(Z(1:100:Nz));
        xlabel('Width (mm)');
        ylabel('Depth (mm)');
        title('Simulation map');

        % Rugosity calculation 
        rugosity = ComputeRugosity(Map, grid, medium);
        format = 'The interface have a rugosity of %.1f %%\n';
        fprintf(format,rugosity)
    end
end

function [pore_x_px, pore_z_px, pore_size_px] = MakePorePorosity(grid, medium, pore, print_hist)
    to_px = @(mm) round(mm/grid.step);

    % PORE AND BONE AREA  
    bone_area = grid.width*medium.interface;
    total_pore_area = bone_area * pore.porosity/100; 
    
    % % PORE SIZE 
    % The pore size follow a gamma distribution, the alpha and beta
    % parameters depends of the patients. 

    if strcmp(pore.patient, 'aged')
        alpha = 25; beta = 2;                       
    elseif strcmp(pore.patient, 'osteoporotic')
        alpha = 10; beta = 7;    
    elseif strcmp(pore.patient, 'young')
        alpha = 35;  beta = 1; 
    else 
        disp('Unknown patient type');
    end
    
    % Generation of diameters until the total area reaches target area
    pore_size = [];                                      % Distribution of the pore diameter (mm)
    pore_area = 0;
    while pore_area <= total_pore_area                   
        pore_size(end+1) = gamrnd(alpha, beta)*1e-3;     % Alpha and beta values corresponds to values in um   
        pore_area = sum(pi * (pore_size / 2).^2);
    end
    
    pore_size_px = to_px(pore_size/2);                   % Distribution of the pore radius (px)
    
    % Plot histogram of simulated diameters
    if print_hist == true
        figure;
        histogram(pore_size);
        xlabel('Diameter (mm)');
        ylabel('Number of pores');
        title('Simulated Diameter Distribution');
    end 
    
    % PORE POSITION
    % Position of the pore along the X-axis 
    pore_nb = length(pore_size);
    pore_x = rand(1, pore_nb)*grid.width;
    pore_x_px = to_px(pore_x);
    
    % Position of the pore is gaussian distributed along the Z-axis,
    % centered at the medium interface
    mu = medium.interface;                          % Center of the distribution (mm)
    std = 0.5*1.1750;                               % FWHM/2 = 0.5 mm for the pore distribution and FWHM = 2.35*std 
    
    pore_z = [];
    % Generate sample
    while length(pore_z) < pore_nb
        sample = normrnd(mu, std);
        % We only want sample in the bone, thus we choose samples under the interface.
        if sample < mu
            pore_z(end+1) = sample;
        end
    end
    
    pore_z_px = to_px(pore_z);
    
    % Plot histogram for verification 
    if print_hist == true
        figure;
        histogram(pore_z);
        xlabel('Position');
        ylabel('Number of pores');
        title('Simulated Position Distribution');
    end
end 

function [rugosity] = ComputeRugosity(Map, grid, medium)
   % The rugosity corresponds to the volume of pore in the 200um width after 
   % bone/soft tissu interface 
   to_px = @(mm) round(mm/grid.step);
   lambda = (medium.cp(2) / (probe.fc*1e3));       % Wavelength in the bone (mm)
   interface = to_px(lambda);
   
   pore_area = sum(Map(to_px(medium.interface) - interface :to_px(medium.interface), :) == 0, 'all');
   tot_area = numel(Map(to_px(medium.interface) - interface :to_px(medium.interface), :));
   
   rugosity = pore_area / tot_area *100;
end