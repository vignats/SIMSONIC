function [Rq, Corr, rugosity] = ComputeInterfaceParameters(Map, grid, probe, medium)
% This function compute the parameters that allows to describe the
% interface. 
% INPUT :   Map - matrice of 0 and 1 that describe the interface between the
% bone and soft tissus.
%           gird - structur that contain the information about the step of
%           the map, in mm.
%           probe - structur that describe that contain the information 
%           about the central frequency of the probe, in Hz. 
%           medium - structur that describe that contain the information 
%           about the velocity in the different mediums, in m/s.
%
% OUTPUT :  Rq - Root mean square of the height profile of the interface.
%           Corr - Correlation lenght of the height profile of the interface.
%           rugosity - Density of the pore in a region of 1 wavelength
%           width around the porous interface.
%
% See also : ComputeCorr
 
    % Extract the bone interface
    profile = zeros(1, size(Map, 2));
    for i = 1:size(Map, 2)
        profile(1, i) = find(Map(:, i) == 1, 1, 'last'); % Find last bone element that correspond to the bone ST interface
    end
    interfacePorous = max(profile);                % In pixel. 

    % Computation of the parameters for a profile centered around zero.
    profile = profile - mean(profile);
    Rq = rms(profile)*grid.step;                    % RMS of the surface (mm)      
    Corr = ComputeCorr(profile, grid.step);         % Correlation length of the surface (mm)

    % The rugosity corresponds to the volume of pore in one wavelength width after 
    % bone/soft tissu interface 
    to_px = @(mm) round(mm/grid.step);
    lambda = (medium.cp(2) / (probe.fc))*1e3;       % Wavelength in the bone (mm)
    lambda = to_px(lambda);
    
    pore_area = sum(Map(round(interfacePorous) - lambda :round(interfacePorous), :) == 0, 'all');
    tot_area = numel(Map(round(interfacePorous) - lambda :round(interfacePorous), :));
    
    rugosity = pore_area / tot_area *100;
   
end
