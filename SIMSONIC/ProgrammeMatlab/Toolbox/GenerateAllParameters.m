function [param, grid, probe, medium, interface, signal, simu_dir] = GenerateAllParameters(saveParam)
% This function contain all the parameters required for the simulation. 
% Needs to be generate in MakeFiles.m
% 
% Input : 
% In order to save in a directory a .mat file that contains all the

    % ACQUISITION PARAMETERS
    param.length = 20;               % Duration of the simulation (us) = sqrt((medium.interface - probe.depth)^2 + grid.width^2)*2e3/max(medium.cp) (us)
    param.cfl = 0.99;                % CFL 
    param.PML = 3;                   % PML thicknes (mm)
    param.PML_eff = 120;             % PML efficiencty (dB)
    param.Vmax = 3500;               % Maximal velocity (m/s) = max(medium.cp)
    
    % GRID PARAMETERS 
    grid.step = 0.01;                % Grid step (mm)
    grid.depth = 15;                 % Grid depth or Z-axis (mm)
    grid.width = 300;                 % Grid with or X-axis (mm)
    
    % PROBE PARAMETERS 
    probe.depth = 2;                 % Probe depth (mm)
    probe.fc = 2.5e6;                % Central frequency (Hz)
    probe.pitch = 0.3;               % Pitch (mm)
    probe.width = 0.25;              % Width of the element (mm)
    probe.Nelements = 96;            % Number of elements 
    
    % MEDIUM PARAMETERS
    medium = struct;
    medium.cp(1) = 1540;             % Longitudinal velocity (in m/s, default = 1540 m/s)
    medium.cp(2) = 3500;     
    medium.cs(1) = 0;                % Shear velocity (in m/s, default = 0 m/s)
    medium.cs(2) = 0;        
    medium.rho(1) = 1000;            % Density (kg/m3)
    medium.rho(2) = 2000;     
    medium.attenuation(1) = 0;       % Attenuation coefficient (dB/cm/MHz, default: 0)
    medium.attenuation(2) = 0;     
    
    % INTERFACE PARAMETERS
    interface.depth = 10;                 % Interface between the bone and the soft tisse (mm)
    interface.patient = 'osteoporotic';   % Type of patient that determine the pore size distribution ('young', 'aged' or 'osteoporotic')
    interface.corr = 0.5;                 % correlation length (mm)
    interface.rms = 0.2;                % rms height (mm)
    interface.porosity = 50;              % Porosity in the bone (volume of pore/volume of bone) (%)
    interface.rugosity = 60;              % Rugosity in a layer of a wavelength size at the bone interface (%)
    
    % SIGNAL PARAMETERS 
    signal.fc = probe.fc;            % Central frequency (Hz)
    signal.B = 1.33e6;               % Bandwith at 3dB (Hz)
    signal.nb_periode = 3;           % Number of periode for the emitted signal

    % FILTRATION PARAMETERS
    % To create interface based on the waviness and/or roughness of the
    % endost of a ex-vivo bone endost imaged by X-Ray, contained in the RMS
    % directory. The waviness and/or roughness are obtained by filtration
    % of the endost boundary in the frequency domain

    % CREATION OF THE SIMULATION DIRECTORY
    simuDirAll = '/calculSSD/salome/Simulation-28mars';
    format = 'simulation_rms_%.2f_cl_%.1f/';
    simulation_name = sprintf(format, interface.rms, interface.corr);
    simu_dir = [simuDirAll, simulation_name]; 
    if ~exist(simu_dir,'dir')
        mkdir(simu_dir);
    end

    % SAVE PARAMETERS AND USEFUL FUNCTION  
    if saveParam
        save(fullfile(simu_dir, 'parameters.mat'))
    end
end