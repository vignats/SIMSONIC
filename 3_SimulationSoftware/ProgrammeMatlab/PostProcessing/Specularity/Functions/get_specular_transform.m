% [SPECULAR_TRANSFORM_BONE_RAW] = get_specular_transform(...
%     RECON_PARAM,TOF_BONE,AOV_elem_deg_BONE, PROBE_PARAM.FS, TILT_angles);
% SPECULAR_TRANSFORM_BONE = SPECULAR_TRANSFORM_BONE_RAW;%SPECULAR_TRANSFORM_BONE_INTERP;%

function [SPECULAR_TRANSFORM] = get_specular_transform(...
    RECON_PARAM,TOF,AOV, FS, TILT_angles)
% TODO: Add function description
    arguments
        RECON_PARAM struct
        TOF struct
        AOV struct % incident and reflected angle at pixel
        FS (1,1) 
        TILT_angles {mustBeInRange(TILT_angles,-90,90)}
    end   
    [~,Nz,Nx] = size(TOF.Time_T);
    NTilts  = numel(TILT_angles);
    SIG = RECON_PARAM.rf_data;
    SIG = SIG(RECON_PARAM.PROBE_PARAM.offset+1:end,:,:);
    % I/Q separation.
    I_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
    Q_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
    for Tx = 1:size(SIG,3)
        tmp =  hilbert(squeeze(SIG(:,:,Tx)));
        I_SIG(:,:,Tx) = real(tmp);
        Q_SIG(:,:,Tx) = imag(tmp);
    end
    IQ_SIG = complex(I_SIG,Q_SIG);
    [Nsamples, NRx, NTx] = size(IQ_SIG);
    Time_T = TOF.Time_T;        Time_R = TOF.Time_R;
    Angle_T=AOV.Angle_T;        Angle_R=AOV.Angle_R;
    SPECULAR_TRANSFORM = zeros(NTilts,Nz,Nx);
    fprintf('-------------\n Specular transformation \n');
    % wb=waitbar(0,'1D specular transformation ...');
    for iz=1:Nz
        % waitbar(iz/Nz,wb,'1D specuar beamforming ...');
        for ix =1:Nx
            tx_angled = Angle_T(:,iz,ix);rx_angled = Angle_R(:,iz,ix);
            tx_tflight = Time_T(:,iz,ix);rx_tflight = Time_R(:,iz,ix);
            for it=1:NTx
                at = tx_angled(it);
                del_t = tx_tflight(it);
                if ~del_t
                    continue
                end                
                for ir=1:NRx
                    ar = rx_angled(ir);
                    delay = (del_t+rx_tflight(ir))*FS+1;
                    idelay = floor(delay);
                    weight = delay-idelay;
                    for ithet=1:NTilts
                        thet = TILT_angles(ithet);
                        if abs(ar+at-2*thet)<=1
%                             linidelay_1=idelay+(ir-1)*Nsamples+...
%                                 (itx-1)*Nsamples*NRx;
%                             disp([iz ix]);
                            if idelay<Nsamples
                                val= (1-weight).*IQ_SIG(idelay,ir,it) ...
                                     + weight.*IQ_SIG(idelay+1,ir,it);
                                SPECULAR_TRANSFORM(ithet,iz,ix) = ...
                                    SPECULAR_TRANSFORM(ithet,iz,ix)+val;

                            end
                        end
                    end
                end
            end

        end
    end
    % close(wb)
end

