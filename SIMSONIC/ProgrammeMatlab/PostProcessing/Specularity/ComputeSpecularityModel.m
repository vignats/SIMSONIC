function [SPECULAR_MODEL] = ComputeSpecularityModel(SPECULAR_TRANSFORM, acquisition, reconstruction, TiltAngles, plot, simu_dir)
    % MODEL OF SPECULAR TRANSFORM : simplified version
    model_version = 'exact_1';%exact_0 exact_1 simplified_0 simplified_1
    excitation_signal=hilbert(SimSonic2DReadSgl([simu_dir '/Signal.sgl']));
    
    acquisition.XR = acquisition.XS;
    fprintf('---------GET SPECULAR TRANSFORM MODELS---------\n')
    
    fprintf('\tSoft Tissues\n')
    tic
    % MODEL_MAP_TISSUE = get_simplified_model_1D_tissu_NEW(estimated_geometry,PROBE_PARAM, ...
    %     TIME_of_FLIGHT,excitation_signal,PROBE_PARAM.FS,ANGLE_of_VIEW,TiltAngles);
    NZ = numel(reconstruction.Z); NX = numel(reconstruction.X);
    MODEL_MAP = zeros([length(TiltAngles) NZ NX]);
    SPECULAR_INTERFACE = [0 0]; % coefficient d'une droite [a b] z=ax+b
    wb=waitbar(0,'Simplified specuar model...');
    for iz=1:NZ
        waitbar(iz/NZ,wb,'Simplified specular model ...');
        zp = reconstruction.Z(iz);
        for ix=1:NX
            transmit_time = reconstruction.timeFlight.Time_T(:,iz,ix);
            transmit_angle = reconstruction.angleView.Angle_R(:,iz,ix);
            xp = reconstruction.X(ix);
            if all(transmit_time)
                SPECULAR_MODEL=function_get_pix_simplified_tissue_model(xp,zp,...
                    SPECULAR_INTERFACE,reconstruction.BEST_C_LAYER,acquisition.XS, ...
                    excitation_signal,acquisition.Fs);
                % Interpolation to match desired specular angles
                SPECULAR_MODEL = interp1(transmit_angle,SPECULAR_MODEL,TiltAngles.');
                MODEL_MAP(:,iz,ix) = SPECULAR_MODEL;
            end
        end
    end
    toc
    close(wb)

    if plot 
        NTilts = length(TiltAngles);
        SPECULAR_TRANSFORM(isnan(SPECULAR_TRANSFORM))=0;
        MODEL_MAP(isnan(MODEL_MAP)) = 0;
        tic
        CORRELATION = zeros([NZ NX 2*NTilts-1]);
        for iz =1:NZ
            for ix = 1:NX
                x = SPECULAR_TRANSFORM(:,iz,ix);
                y = MODEL_MAP(:,iz,ix);
                CORRELATION(iz,ix,:)  = xcorr(x,y,'normalized');
            end
        end
        toc
        %
        
        thresh = .3;
        D_THETA = (TiltAngles(2)-TiltAngles(1));
        lags = D_THETA*(-NTilts+1:NTilts-1);
        [PROBA_MAP_TISS, ind_T] = max(abs(CORRELATION),[],3);
        PROBA_MAP_TISS(isnan(PROBA_MAP_TISS))=0;
        
        
        SPECULAR_TILT_MAP_TISS = lags(ind_T);
        PROBA_MAP_TISS_bin  = double((PROBA_MAP_TISS>thresh));
        PROBA_MAP_TISS_bin(PROBA_MAP_TISS_bin==0)=nan;
        TILT_MAP_TISS_bin = SPECULAR_TILT_MAP_TISS.*PROBA_MAP_TISS_bin;
        SPECULAR_TILT_MAP_TISS(isnan(PROBA_MAP_TISS_bin))=nan;
        %
        ST = SPECULAR_TILT_MAP_TISS;
        ST(isnan(SPECULAR_TILT_MAP_TISS)) = 0;
        
        figure('Position', [1 1 1244 1321]),
        Xinit = 500;
        Zinit = 600;
        tlo = tiledlayout(1,2);
        % nexttile(tlo,1)
        pcolor(reconstruction.Xmm,reconstruction.Zmm,PROBA_MAP_TISS)
        xlabel('Lateral position [mm]', Interpreter='latex')
        ylabel('Depth [mm]', Interpreter='latex')
        colorbar
        title('Probability map TISSUE',interpreter='latex')
        axis ij image, shading flat
        clim([0 1])
        % hold on
        % plot((Xinit:3000-Xinit)*parameters.grid.step,profile.heights(Xinit:end-Xinit))
        % hold off
        figure('Position',[1317 1 1244 1321] ),
    
        % nexttile(tlo,2)
        pcolor(reconstruction.Xmm, reconstruction.Zmm, SPECULAR_TILT_MAP_TISS) 
        axis tight
        title('Simulation map');   
        shading flat
        xlabel('Width (mm)');
        ylabel('Depth (mm)');
        axis image;
        colorbar
    end
end
