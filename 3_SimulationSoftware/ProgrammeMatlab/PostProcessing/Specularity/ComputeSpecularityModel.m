function [SPECULAR_MODEL, PROBA_MAP, SPECULAR_TILT_MAP] = ComputeSpecularityModel(...
    SPECULAR_TRANSFORM, acquisition, reconstruction, TiltAngles, plot, parameters, simu_dir)
    % MODEL OF SPECULAR TRANSFORM : simplified version
    model_version = 'exact_1';%exact_0 exact_1 simplified_0 simplified_1
    excitation_signal=hilbert(SimSonic2DReadSgl([simu_dir '/Signal.sgl']));
    
    acquisition.XR = acquisition.XS;
    fprintf('---------GET SPECULAR TRANSFORM MODELS---------\n')
    
    if isfield(parameters, 'bone') || (isfield(parameters.interface, 'periost') && parameters.interface.periost > 0)
        % Compute the specular transform in the 3 layers case
        estimated_geometry = struct('X',reconstruction.X,'Z',reconstruction.Z,'lens_thick',acquisition.LENS_THICKNESS,...
         ...
            'fit_curve_Periosteum',reconstruction.Z(reconstruction.PERI(4:end)),...
            'peribola_coef',reconstruction.PERI(1:3),'endobola_coef', reconstruction.ENDO(1:3),...
            'V0',acquisition.C_LENS, ...
            'V1',reconstruction.C_TISSUE, ...
            'V2',reconstruction.C_RADIAL...
            );
        fprintf('\tSoft Tissues\n')
        tic
        SPECULAR_MODEL.Tissu = get_simplified_model_1D_tissu_NEW(estimated_geometry,acquisition, ...
            reconstruction.Tissu.timeFlight,excitation_signal,acquisition.Fs,reconstruction.Tissu.Pixel.angleView,TiltAngles);
        toc

        SPECULAR_MODEL.Tissu(isnan(SPECULAR_MODEL.Tissu)) = 0;
        
        fprintf('\tBone Tissues - 1D model\n')
        tic,
        SPECULAR_MODEL.Bone= get_simplified_model_1D_NEW(estimated_geometry,acquisition, ...
            reconstruction.Bone.timeFlight,excitation_signal,acquisition.Fs,reconstruction.Bone.Pixel.angleView,TiltAngles);
        toc
        SPECULAR_MODEL.Bone(isnan(SPECULAR_MODEL.Bone)) = 0;

        % Compute the correlation
        tic
        CORR3_BONE_nested= xcorr_nested(SPECULAR_TRANSFORM.Bone, ...
            SPECULAR_MODEL.Bone);
        toc
        tic
        CORR3_TISS_nested = xcorr_nested(SPECULAR_TRANSFORM.Tissu, ...
            SPECULAR_MODEL.Tissu);
        toc
        %
        %
        NTilts = length(TiltAngles);
        thresh = .5;
        D_THETA = (TiltAngles(2)-TiltAngles(1));
        lags = D_THETA*(-NTilts+1:NTilts-1);
        [PROBA_MAP.Tissu, ind_T] = max(abs(CORR3_TISS_nested),[],3);
        PROBA_MAP.Tissu(isnan(PROBA_MAP.Tissu))=0;
        [PROBA_MAP.Bone, ind_B] = max(abs(CORR3_BONE_nested),[],3);
        PROBA_MAP.Bone(isnan(PROBA_MAP.Bone))=0;
        
        SPECULAR_TILT_MAP.Tissu = lags(ind_T);
        PROBA_MAP_Tissu_bin  = double((PROBA_MAP.Tissu>thresh));
        PROBA_MAP_Tissu_bin(PROBA_MAP_Tissu_bin==0)=nan;
        TILT_MAP_TISS_bin = SPECULAR_TILT_MAP.Tissu.*PROBA_MAP_Tissu_bin;
        
        SPECULAR_TILT_MAP.Bone = lags(ind_B);
        PROBA_MAP_Bone_bin  = double((PROBA_MAP.Bone>thresh));
        PROBA_MAP_Bone_bin(PROBA_MAP_Bone_bin==0)=nan;
        TILT_MAP_BONE_bin = SPECULAR_TILT_MAP.Bone.*PROBA_MAP_Bone_bin;

    else
        % Compute the specular transform in the 2 layers case
        fprintf('\tBone Tissues\n')
        tic

        NZ = numel(reconstruction.Z); NX = numel(reconstruction.X);
        MODEL_MAP.Bone = zeros([length(TiltAngles) NZ NX]);
        SPECULAR_INTERFACE = [0 0]; % coefficient d'une droite [a b] z=ax+b
        wb=waitbar(0,'Simplified specular model...');
        for iz=1:NZ
            waitbar(iz/NZ,wb,'Simplified specular model ...');
            zp = reconstruction.Z(iz);
            for ix=1:NX
                transmit_time = reconstruction.Bone.timeFlight.Time_T(:,iz,ix);
                transmit_angle = reconstruction.Bone.Degre.angleView.Angle_R(:,iz,ix);
                xp = reconstruction.X(ix);
                if all(transmit_time)
                    SPECULAR_MODEL.Bone=function_get_pix_simplified_tissue_model(xp,zp,...
                        SPECULAR_INTERFACE,reconstruction.C_TISSUE,acquisition.XS, ...
                        excitation_signal,acquisition.Fs);
                    % Interpolation to match desired specular angles
                    SPECULAR_MODEL.Bone = interp1(transmit_angle,SPECULAR_MODEL.Bone,TiltAngles.');
                    MODEL_MAP.Bone(:,iz,ix) = SPECULAR_MODEL.Bone;
                end
            end
        end
        toc
        close(wb)

        % Compute the correlation 
        NTilts = length(TiltAngles);
        SPECULAR_TRANSFORM.Bone(isnan(SPECULAR_TRANSFORM.Bone))=0;
        MODEL_MAP.Bone(isnan(MODEL_MAP.Bone)) = 0;
        tic
        CORRELATION = zeros([NZ NX 2*NTilts-1]);
        for iz =1:NZ
            for ix = 1:NX
                x = SPECULAR_TRANSFORM.Bone(:,iz,ix);
                y = MODEL_MAP.Bone(:,iz,ix);
                CORRELATION(iz,ix,:)  = xcorr(x,y,'normalized');
            end
        end
        toc

        thresh = .3;
        D_THETA = (TiltAngles(2)-TiltAngles(1));
        lags = D_THETA*(-NTilts+1:NTilts-1);
        [PROBA_MAP.Bone, ind_B] = max(abs(CORRELATION),[],3);
        PROBA_MAP.Bone(isnan(PROBA_MAP.Bone))=0;


        SPECULAR_TILT_MAP.Bone = lags(ind_B);
        PROBA_MAP_Bone_bin  = double((PROBA_MAP.Bone>thresh));
        PROBA_MAP_Bone_bin(PROBA_MAP_Bone_bin==0)=nan;
        TILT_MAP_Bone_bin = SPECULAR_TILT_MAP.Bone.*PROBA_MAP_Bone_bin;
        SPECULAR_TILT_MAP.Bone(isnan(PROBA_MAP_Bone_bin))=nan;

    end

    if plot 
        figure('Position', [1 1 1244 1321]),
        Xinit = 500;
        Zinit = 600;
        tlo = tiledlayout(1,2);
        % nexttile(tlo,1)
        pcolor(reconstruction.Xmm,reconstruction.Zmm,PROBA_MAP.Bone)
        xlabel('Lateral position [mm]', Interpreter='latex')
        ylabel('Depth [mm]', Interpreter='latex')
        colorbar
        colormap jet
        title('Probability map',interpreter='latex')
        axis ij image, shading flat
        clim([0 1])
        % hold on
        % plot((Xinit:3000-Xinit)*parameters.grid.step,profile.heights(Xinit:end-Xinit))
        % hold off
        figure('Position',[1317 1 1244 1321] ),
    
        % nexttile(tlo,2)
        pcolor(reconstruction.Xmm, reconstruction.Zmm, SPECULAR_TILT_MAP.Bone) 
        axis tight
        title('Detected specular angle map');   
        shading flat
        xlabel('Width (mm)');
        ylabel('Depth (mm)');
        axis ij image;
        clim([-1 1]*45)
        colorbar
        colormap jet
    end
end
