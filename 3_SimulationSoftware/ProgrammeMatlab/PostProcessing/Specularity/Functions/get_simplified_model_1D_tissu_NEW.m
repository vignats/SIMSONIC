function MODEL_MAP = get_simplified_model_1D_tissu_NEW(estimated_geometry,PROBE_PARAM,...
    TOF,excitation_signal,FS,AOV, TILT_ANGLES)
    XS = PROBE_PARAM.XS;
%     NTX = numel(XS);
    NZ = numel(estimated_geometry.Z);
    NX = numel(estimated_geometry.X);
    NTILTs = numel(TILT_ANGLES);
    MODEL_MAP = zeros([NTILTs NZ NX]);
    SPECULAR_INTERFACE = [0 0];
    wb=waitbar(0,'Simplified specular model...');
    for iz=1:NZ
        waitbar(iz/NZ,wb,'Simplified specular model ...');
        for ix=1:NX
            transmit_time = TOF.Time_T(:,iz,ix);
            transmit_angle = AOV.Angle_R(:,iz,ix);
            if all(transmit_time)
%                 SPECU_MODEL_MATLAB_BONEE = zeros([NTILTs 1]);
                SPECU_MODEL_MATLAB_BONEE=get_pix_simplified_tissue_model_NEW(ix,iz,...
                    SPECULAR_INTERFACE,estimated_geometry,XS,excitation_signal,transmit_time,FS);
                SPECU_MODEL_MATLAB_BONEE = interp1(transmit_angle,SPECU_MODEL_MATLAB_BONEE,TILT_ANGLES.');
                MODEL_MAP(:,iz,ix) = SPECU_MODEL_MATLAB_BONEE;
            end
        end
    end
    close(wb);
end