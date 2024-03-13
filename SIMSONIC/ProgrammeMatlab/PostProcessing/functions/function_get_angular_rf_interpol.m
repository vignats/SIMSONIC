function [RF_MAP] = function_get_angular_rf_interpol(...
    RECON_PARAM,TOF,FS)
% TODO: Add function description
    arguments
        RECON_PARAM struct
        TOF (1,1) struct
        FS (1,1) 
    end   
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
    Time_T = TOF.Time_T;Time_R = TOF.Time_R;
    [~,Nz,Nx] = size(Time_T);
    offset_lin_idx=(0:NRx-1)*Nsamples+...
                (0:NTx-1)'*Nsamples*NRx;
    RF_MAP = nan(NRx,NTx,Nz,Nx,'double');%+0*1i;
    wb=waitbar(0,'Getting map of rf data ...');
    for iz=1:Nz%
        p = iz/Nz;
        waitbar(p,wb,['Getting map of rf data ... ' ...
            num2str(round(p*100)) '%']);
        if ~nnz(Time_T(:,iz,:))
            continue
        end
        for ix=1:Nx%
            tx_tflight = Time_T(:,iz,ix);rx_tflight = Time_R(:,iz,ix);
            if ~nnz(tx_tflight)
                continue
            end            
            delay_sample = (tx_tflight'+rx_tflight)*FS+1; %[NRx x NTx]
            clear tx_tflight rx_tflight
            idelay_sample= floor(delay_sample); 
            idelay_sample(idelay_sample>=Nsamples)=nan;
%             m=find(~isnan(idelay_sample));
            lin_idelay=idelay_sample+offset_lin_idx; % Linear index
            weight= delay_sample-idelay_sample;
%             weight(isnan(weight))=[];
            lin_idelay(isnan(lin_idelay)) = 1;%[];
            clear idelay_sample delay_sample
            RF_MAP(:,:,iz,ix) =(1-weight).*IQ_SIG(lin_idelay) ...
                + weight.*IQ_SIG(lin_idelay+1);
%             RF_MAP(sub2ind([NRx*NTx,Nz,Nx],m,iz*ones(size(m)),ix*ones(size(m)))) = (1-weight).*IQ_SIG(lin_idelay) ...
%                             + weight.*IQ_SIG(lin_idelay+1);
            clear weight lin_idelay 
        end
    end
    close(wb)
end

