function [param, grid, probe, medium, signal] = GenerateAllParametersExVivo()
% This function contain all the parameters required for the simulation. 
% Needs to be generate in MakeFiles.m
% 
% INPUT : saveParam - boolean that indicate weither to save the parameters
% in a .mat file in the simulation directory.
  
    % FILTRATION PARAMETERS
    % To create interface based on the waviness and/or roughness of the
    % endost of a ex-vivo bone endost imaged by X-Ray, contained in the RMS
    % directory. The waviness and/or roughness are obtained by filtration
    % of the endost boundary in the frequency domain
    
    % ACQUISITION PARAMETERS
    param.length = 20;               % Duration of the simulation (us) = sqrt((medium.interface - probe.depth)^2 + grid.width^2)*2e3/max(medium.cp) (us)
    param.cfl = 0.99;                % CFL 
    param.PML = 3;                   % PML thicknes (mm)
    param.PML_eff = 120;             % PML efficiencty (dB)
    param.Vmax = 3500;               % Maximal velo1.25city (m/s) = max(medium.cp)
    
    % GRID PARAMETERS 
    grid.step = 0.01;                % Grid step (mm)
    grid.depth = 15;                 % Grid depth or Z-axis (mm)
    grid.width = 30;                 % Grid with or X-axis (mm)
    
    % PROBE PARAMETERS 
    probe.depth = 2;                 % Probe depth (mm)
    probe.fc = 2.5e6;                % Central frequency (Hz)
    probe.pitch = 0.3;               % Pitch (mm)
    probe.width = 0.25;              % Width of the element (mm)
    probe.Nelements = 96;            % Number of elements 
    
    % MEDIUM PARAMETERS
    medium = struct;    bone.fc = 0.0714;         % Cut-off frequency, fc = 0.06 for waviness + roughness and 1.25 for roughness (mm-1)
    bone.fs = 9e-3;           % Pixel size of the X-Ray image (mm)
    bone.segmented = false;   % Indicate if the image is segmented
    medium.cp(1) = 1540;             % Longitudinal velocity (in m/s, default = 1540 m/s)
    medium.cp(2) = 3500;     
    medium.cs(1) = 0;                % Shear velocity (in m/s, default = 0 m/s)
    medium.cs(2) = 0;        
    medium.rho(1) = 1000;            % Density (kg/m3)
    medium.rho(2) = 2000;     
    medium.attenuation(1) = 0;       % Attenuation coefficient (dB/cm/MHz, default: 0)
    medium.attenuation(2) = 0;     
      
    % SIGNAL PARAMETERS 
    signal.fc = probe.fc;            % Central frequency (Hz)
    signal.B = 1.33e6;               % Bandwith at 3dB (Hz)
    signal.nb_periode = 3;           % Number of periode for the emitted signal

end