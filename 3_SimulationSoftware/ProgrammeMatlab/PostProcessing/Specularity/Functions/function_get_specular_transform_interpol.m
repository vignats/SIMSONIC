function [SPECULAR_TRANSFORM] = function_get_specular_transform_interpol(...
    RECON_PARAM,TOF,AOV_px,AOV_elem, FS, TILT_angles)
% TODO: Add function description
    arguments
        RECON_PARAM struct
        TOF struct
        AOV_px struct % incident and reflected angle at pixel
        AOV_elem struct%
        FS (1,1) 
        TILT_angles {mustBeInRange(TILT_angles,-90,90)}
    end
% TODO: Add function description
    % TOF = RECON_PARAM.timeFlight;
    % AOV_px = RECON_PARAM.angleView;
    % AOV_elem = AOV_px;
    % FS = acquisition.Fs;

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
%     lambda_eau = 1540/2.5e6;
%     directivity = @(theta) sinc(295e-6/lambda_eau*sind(theta)).*cosd(theta);
    max_half_open_angle = 180;
    [Nsamples, NRx, NTx] = size(IQ_SIG);
    [~,Nz,Nx] = size(TOF.Time_T);
    NTilts  = numel(TILT_angles);
    Time_T = TOF.Time_T;        Time_R = TOF.Time_R;
    Angle_T_px=AOV_px.Angle_T;        Angle_R_px=AOV_px.Angle_R;
    Angle_T_elem=AOV_elem.Angle_T;        Angle_R_elem=AOV_elem.Angle_R;
    SPECULAR_TRANSFORM = zeros(NTilts,Nz,Nx);
    wb=waitbar(0,'1D specular transformation ...');
    for iz=1:Nz%
        waitbar(iz/Nz,wb,'1D specular transformation ...');
        for ix=1:Nx%
            tx_tflight = Time_T(:,iz,ix);rx_tflight = Time_R(:,iz,ix)';
            tx_angled = Angle_T_px(:,iz,ix);rx_angled = Angle_R_px(:,iz,ix); 
            tx_angled_elem = Angle_T_elem(:,iz,ix);rx_angled_elem = Angle_R_elem(:,iz,ix); 
            tmp = SPECULAR_TRANSFORM(:,iz,ix);
            for itx=1:NTx
                emit_directivity = 1;%directivity(tx_angled_elem(itx));
                % For each tilt angle, get corresponding specular angle
                specular_angle = -tx_angled(itx)-2*TILT_angles; % 
                del_t = tx_tflight(itx);
                % Pour chaque tilt trouver l'indice du 
                % récepteur qui minimise
                % |alpha_r+alpha_t-2Theta|>0
                if del_t==0 || abs(tx_angled_elem(itx))>max_half_open_angle
                    continue
                end
                idx_nearest_rx_1 = get_idx_of_closest_value_v2(rx_angled,specular_angle);
                idx_nearest_rx_2 = idx_nearest_rx_1+...
                    sign(rx_angled(idx_nearest_rx_1)'-specular_angle);%*1;
%                 idx_nearest_rx_2(idx_nearest_rx_2>NTx)=NTx;
%                 idx_nearest_rx_2(idx_nearest_rx_2<1)=1;
                good_tilt_idx = (idx_nearest_rx_2<=NTx & idx_nearest_rx_2>=1);

%                 good_tilt_idx = (idx_nearest_rx_1~=idx_nearest_rx_2);
                relevant_specular_angle=specular_angle(good_tilt_idx)';                
                idx_nearest_rx_1 = idx_nearest_rx_1(good_tilt_idx);
                idx_nearest_rx_2 = idx_nearest_rx_2(good_tilt_idx);
            
                nearest_rxangle1 = rx_angled(idx_nearest_rx_1); % NTilts X 1
                nearest_rxangle2 = rx_angled(idx_nearest_rx_2);
                delay_1= (del_t+rx_tflight(idx_nearest_rx_1))*FS+1;
                delay_2= (del_t+rx_tflight(idx_nearest_rx_2))*FS+1;
                idelay_1= floor(delay_1); idelay_2= floor(delay_2);

                weight_1= delay_1-idelay_1; weight_2= delay_2-idelay_2;
                linidelay_1=idelay_1+(idx_nearest_rx_1-1)*Nsamples+...
                    (itx-1)*Nsamples*NRx;
                linidelay_2=idelay_2+(idx_nearest_rx_2-1)*Nsamples+...
                    (itx-1)*Nsamples*NRx;
                poss_1 = idelay_1<Nsamples;
                poss_2 = idelay_2<Nsamples;
                rx_angle_elem_1 = rx_angled_elem(idx_nearest_rx_1(poss_1));
                rx_angle_elem_2 = rx_angled_elem(idx_nearest_rx_2(poss_2));
                receive_directivity_1 = ones(size(rx_angle_elem_1'));%directivity(rx_angled_elem(idx_nearest_rx_1(poss_1))).';
                receive_directivity_2 = ones(size(rx_angle_elem_2'));%directivity(rx_angled_elem(idx_nearest_rx_2(poss_2))).';                    
                receive_directivity_1(abs(rx_angle_elem_1)>max_half_open_angle)=0;
                receive_directivity_2(abs(rx_angle_elem_2)>max_half_open_angle)=0;
                linidelay_1(~poss_1) = []; linidelay_2(~poss_2) = [];
                weight_1(~poss_1) = []; weight_2(~poss_2) = [];
                val_1 = zeros(size(idelay_1)); val_2 = zeros(size(idelay_2));
                val_1(poss_1)= emit_directivity*(1-weight_1).*receive_directivity_1.*IQ_SIG(linidelay_1) ...
                    + emit_directivity*weight_1.*receive_directivity_1.*IQ_SIG(linidelay_1+1); %  [1 x NTilts ]
                val_2(poss_2)= (1-weight_2).*receive_directivity_2.*IQ_SIG(linidelay_2)...
                    + weight_2.*receive_directivity_2.*IQ_SIG(linidelay_2+1);
                weight = (relevant_specular_angle-nearest_rxangle1)...
                    ./(nearest_rxangle2-nearest_rxangle1);% % NTilts X 1
%                 weight(weight>1 | weight<0) = Inf;
                weight = weight.';
                val_interp = (1-weight).*val_1+weight.*val_2; %1 x NTilts 
                val_interp(weight>1 | weight<0)=0;
                tmp(good_tilt_idx) = ...
                    tmp(good_tilt_idx) +val_interp.';
            end
            SPECULAR_TRANSFORM(:,iz,ix) = tmp;
        end
    end
    close (wb)
end

