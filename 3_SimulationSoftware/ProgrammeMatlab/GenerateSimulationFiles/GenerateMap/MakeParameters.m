function[] = MakeParameters(param, grid, probe, medium, simu_dir, print)
% MakeParameters generates the Parameters.ini2D file for the simulation in SimSonic, for each transmitter, reiceiver.
% It creates a directory for each transmitter, and copy the signal and
% geometry file in it. 
% The functions need to be used in the MakeFiles.m script
% 
% Syntax:
%   MakeParameters(param, grid, probe, medium, simu_dir, print)

    general_param = WriteGeneralProperties(param, grid, print);
    bound_pml = WriteBoundaryConditions(param, grid, print);
    snapshot = WriteSnapshotProperties(print);
    medium_str = WriteMaterialProperties(probe,medium, print);

    for tx = 1:probe.Nelements
        % Computation of the parameters for each transmitters.
        % Some parameters about the emission (delay, apodisation, number of
        % elements per emitters...) can be modified in the function 
        % WriteEmitterProperties() of the MakeParameters script. 
        % tx correspond to the emitter number at each position
    
        % Writting of the Parameters.init2D file for each transmitter
        [receivers, emitters] = WriteSensorProperties(tx, grid, probe, print);
    
        all = [general_param newline bound_pml newline ...
                 '%%%%%%% SENSORS %%%%%%%' emitters.T11' emitters.T22'...
                 receivers.T11' ...receivers.T22' receivers.T12'...
                 receivers.V1' newline ...receivers.V2' newline...
                 snapshot newline...
                 medium_str...
                 ];
    
        % Print in each directory
        tx_dir = sprintf('%stx_%02d/',simu_dir, tx);
        if ~exist(tx_dir,'dir')
		    mkdir(tx_dir);
        end
        output_param_filename = [tx_dir 'Parameters.ini2D'];
        fid = fopen(output_param_filename, 'w');
    
        for i = 1:numel(all)
            if i == numel(all)
                fprintf(fid,'%s', all{i});
                break
            else
                fprintf(fid,'%s\n', all{i});
            end
        end
        fclose(fid);
    
        % Copy the signal and Geometry file in each 
        copyfile([simu_dir, '/Signal.sgl'], [tx_dir, 'Signal.sgl']);
        copyfile([simu_dir, '/Geometry.map2D'], [tx_dir, '/Geometry.map2D']);
    
    end

end 

function[general_param] = WriteGeneralProperties(param, grid, print)
    % Writting
    grid_step_str = sprintf('%-29s: %-2.2g','Grid Step           (mm)',grid.step);
    CFL_coef = sprintf('%-29s: %-2.2g','CFL Coefficient',param.cfl);
    VMAX = sprintf('%-29s: %-2.4f','Vmax            (mm/us)',param.Vmax*1e-3);
    simul_length = sprintf('%-29s: %-d','Simulation Length   (us)',round(param.length));
    absorption_type=sprintf('%-29s: %-d\n','Absorption Type',2);%'';

    general_param = {'%%%%%%%%%% GENERAL PARAMETERS %%%%%%%' ...
     grid_step_str' VMAX' CFL_coef' simul_length' absorption_type'};
    
    % Print for verification
    if print == true
        for i=1:length(general_param)
            fprintf('%s\n', general_param{i})
        end
    end
end

function[bound_pml] = WriteBoundaryConditions(param, grid, print)
    % Parameters
    PML_th_grid = param.PML/grid.step;
    PML_eff = param.PML_eff;
    
    % Writting
    boundaries=sprintf('%-29s: %-d\n','X1_low Boundary',0,...
                                      'X1_high Boundary',0,...
                                      'X2_low Boundary',0,...
                                      'X2_high Boundary',0);
    
    PML_thicknes = sprintf('%-29s: %-4g','PML Thickness   (grid step)',PML_th_grid);
    PML_eff = sprintf('%-29s: %-.1f','PML Efficiency     (dB)',PML_eff);
    Vmax_pml = sprintf('%-29s: %-2.4f','Vmax in PML     (mm/us)',param.Vmax*1e-3);
    bound_pml = {'%%%%%%% BOUNDARY COUNDITIONS %%%%%%%' ...
        boundaries PML_thicknes Vmax_pml PML_eff};
    
    % Print for verification
    if print == true
        for i=1:length(bound_pml)
            fprintf('%s\n', bound_pml{i})
        end
    end
end

function[receivers, emitters] = WriteSensorProperties(tx, grid, probe, print)
    to_px = @(mm) round(mm/grid.step);              % To change coordinate in mm in pixel position.
    
    % Sensor parameters (receivers + emitters)
    elt_width_px = to_px(probe.width);    % 240 um to mimic the P4-1 probe 245 um
    picth_px = to_px(probe.pitch);       

    X = 0:grid.step:grid.width-grid.step;       % Position (mm)
    lateral_center = mean(X);                   % Center of the X-axis 
    
    % Sensors position
    position_in_probe = 0:probe.pitch:(probe.Nelements-1)*probe.pitch;   % Position of each element in the probe (mm)
    position_in_probe = position_in_probe - mean(position_in_probe);     % Position of each element regarding the center of the probe (mm)
    position_in_grid = position_in_probe + lateral_center;               % Position of the element from the left side of the grid (mm)
    assert(all(position_in_grid>0));

    sensors_disposition_px = to_px(position_in_grid);                    % Position of the element from the left side of the grid (pixel) 

    elt_z = round(probe.depth/grid.step);                                % Position of the array on the Z-axis (pixel)
    elt_x0 = sensors_disposition_px(1);                                  % Position of the array on the X-axis (pixel)

    % Writting for receivers
    rcv_format_str = '%-29s: %d\n%s\n1\n%-6d%-d\n%-6d%-6d%-d';
    receivers.V1=sprintf(rcv_format_str,'Nb of V1 Receivers Arrays',...
        1,'V1_main',elt_z,elt_x0,probe.Nelements,picth_px,elt_width_px);
    receivers.V2=sprintf(rcv_format_str,'Nb of V2 Receivers Arrays',...
        1,'V2_main',elt_z,elt_x0,probe.Nelements,picth_px,elt_width_px);
    receivers.T11=sprintf(rcv_format_str,'Nb of T11 Receivers Arrays',...
        1,'T11_main',elt_z,elt_x0,probe.Nelements,picth_px,elt_width_px);
    receivers.T22=sprintf(rcv_format_str,'Nb of T22 Receivers Arrays',...
        1,'T22_main',elt_z,elt_x0,probe.Nelements,picth_px,elt_width_px);
    receivers.T12=sprintf(rcv_format_str,'Nb of T12 Receivers Arrays',...
        1,'T12_main',elt_z,elt_x0,probe.Nelements,picth_px,elt_width_px);

    % Print for verification
    if print == true
        fprintf('%s\n', receivers.T11);
        fprintf('%s\n', receivers.V1);
    end
    
    % Transmitter parameters
    nb_tx_per_emission = 1;         % Number of transduceur per emission      
    apo=0;                          % Apodisation, 0 if one element, 1 for Hann window apodisation
    speed_del = 0;                  % Veloxity to compute the delay for focusing or deflecting (mm/us)
    deflection = 0;                 % Emission angle for plane wave (degrees)
    pitch_tx = 0;
    
    tx_x = uint16(sensors_disposition_px(tx));  % Position of the tx sensor (pixel)

    % Writting
    tx_format_str = '%-29s: %d\n-1  Signal.sgl\n1\n%-6d%-d\n%-6d%-6d%-6d%-6d%-6d\n%-2.1d%5.1d';
    emitters.T11=sprintf(tx_format_str,'Nb of T11 Emitters Arrays',...
            1,...
            elt_z,tx_x,...
            nb_tx_per_emission,pitch_tx,elt_width_px,apo,0,...
            deflection,speed_del);
    emitters.T22=sprintf(tx_format_str,'Nb of T22 Emitters Arrays',...
            1,...
            elt_z,tx_x,...
            nb_tx_per_emission,pitch_tx,elt_width_px,apo,0,...
            deflection,speed_del); 

    % Print for verification
    if print == true
        fprintf('%s\n', emitters.T11);
    end
end 

function[snapshot] = WriteSnapshotProperties(print)
    % Writting
    snapshot_rec = sprintf('Snapshots Record Period  (us): .5');
    snap_V1  = sprintf('%-29s: %-d','Record V1 Snapshots',0);
    snap_V2  = sprintf('%-29s: %-d','Record V2 Snapshots',0);
    snap_V   = sprintf('%-29s: %-d','Record V Snapshots',0);
    snap_T11 = sprintf('%-29s: %-d','Record T11 Snapshots',0);
    snap_T12 = sprintf('%-29s: %-d','Record T12 Snapshots',0);
    snapshot = {'%%%%%%%%%% SNAPSHOT PARAMETERS %%%%%%%%%% ' ...
        snapshot_rec' snap_V' snap_T11' };
    
    % Print for verification
    if print == true
        for i=1:length(snapshot)
        fprintf('%s\n', snapshot{i});
        end
    end
end 

function[medium_str] = WriteMaterialProperties(probe, medium, print)
    % Format 
    media{1} = sprintf('%-8s%-8s%-8s%-8s%-8s%-7s%-4s%-4s%-4s%-4s%-10s%-2s', ...
    'Indix','Density', 'c11','c22','c12','c66',...
    'x11','x22','x12','x66','Qval','speed');
    medium_format_str='%-8d%-8.1f%-8.3g%-8.3g%-8.3g%-7.1f%-4d%-4d%-4d%-4d%-10.4g%-2.5g';

    % Parameters    
    k=1;
    x11=0;x22=0;x12=0;x66=0;
    
    for i=1:length(medium.cp)
        k=k+1;
        % Scaled velocities (in mm/us)
        cp = medium.cp(i)*1e-3;
        cs = medium.cs(i)*1e-3;
        % Scaled density (in g/m3)
        rho = medium.rho(i)*1e-3;
        % Elasticity coefficient computation
        c11 = rho*cp^2;
        c22 = c11;
        c66 = rho*cs^2;
        c12 = c11 - 2*c66;

        media{k}=sprintf(medium_format_str,i-1,rho,...
        c11,c22,c12,c66,...
        x11,x22,x12,x66,medium.attenuation(i)*probe.fc*1e-6/10,cp);   
    end

    medium_str = {'%%%%%%%%%% DEFINITION OF MATERIALS PROPERTIES %%%%%%%%%%' ...
    media{1} newline ...
    'Starts Materials List' media{2:end} 'Ends Materials List'};

    % Print for verification
    if print == true
        for elt=medium_str
            fprintf('%s\n',elt{:})
        end
    end
end
