function [BoneSpeed] = Autofocus(RECON_PARAM, PROBE_PARAM, plot1, plot2)
%--- IMPORTANT ---%
% ALL PARAMETERS MUST BE OF CLASS REAL DOUBLE
% VECTORS MUST BE ROW VECTORS
%
%--- WARNING ---%
% the autofocus approach requires a rather larger number of different
% transmit beams (different insonofication angles), typically at least 20
%
%---- Input parameters -----%
% LENS_THICKNESS: Lens thickness (meter)
% C_LENS: known wavespeed in the lens of the probe (m/s)
% FS: Sampling frequency (Hz)
% HalfOpeningAngInLensRad: max half opening angle at receive element (radian)
% C_TISSUE: a priori wavespeed in soft tissue between probe and bone (m/s)
% NCORE: number of CPU cores to use for OMP parallel computing
% (X_El, Z_El): coordinates of array elements (meter)
% (XS, ZS): coordinates of (virtual) sources (meter)
% add_to_delay_firing: travel time from virtual sources to nearest array element (second), zeros for single-element transmissions
% (X, Z): coordinates of image pixels (meter)
% (I_SIG, Q_SIG): in-phase and quadrature sampled RF data, 3D œarray dimension is (Time,Receiver index,Transmission index)
% Tx_to_be_used: indices of transmissions to be used
% C_TISSUE_TEST_VALUES: row vector containing values you want to test

    % Segmentation
    C_TISSUE = RECON_PARAM.C_TISSUE;
    fine_resolution = C_TISSUE/(PROBE_PARAM.FREQ_Transducteur*12); % [m]
    
    XMinSegmentationEndo = PROBE_PARAM.xmin;
    XMaxSegmentationEndo = PROBE_PARAM.xmax;
    %

    % create image axes/pixels coordinates
    X       = PROBE_PARAM.xmin:fine_resolution:PROBE_PARAM.xmax; % image width [m]
    Z       = PROBE_PARAM.zmin:fine_resolution:PROBE_PARAM.zmax; % image depth [m]
    
    C_AXIAL = C_TISSUE;
    C_RADIAL = C_TISSUE;
  
    Setup = [PROBE_PARAM.LENS_THICKNESS RECON_PARAM.D_PROBE_PERIOS...
        RECON_PARAM.MIN_CORTICAL RECON_PARAM.MAX_CORTICAL ...
        RECON_PARAM.HalfOpeningAngInLensRad RECON_PARAM.HalfOpeningAngInLensRad...
        RECON_PARAM.FnumberMin RECON_PARAM.Max_err_allowed_ReceiveAngle_rad PROBE_PARAM.Fs ...
        PROBE_PARAM.C_LENS C_TISSUE RECON_PARAM.C_MARROW...
        C_AXIAL C_RADIAL RECON_PARAM.ANISO_SHAPE_COEF RECON_PARAM.NCORE...
        RECON_PARAM.ReConTo RECON_PARAM.NeedTravelTime RECON_PARAM.SubApertureApodis...
        XMinSegmentationEndo XMaxSegmentationEndo];

    %   
    displayed_dynamic_range_dB = 40;
    PeriParabIn = [0 0 0];EndoParabIn=[0 0 0];
    ReceiveAngAtPix_rad = deg2rad([0]);
    SIG = RECON_PARAM.rf_data(PROBE_PARAM.offset+1:end,:,:);
    Tx_to_be_used = 1:PROBE_PARAM.NSOURCES;
    
    add_to_delay_firing = zeros(1,PROBE_PARAM.NSOURCES);
    % I/Q separation.
    I_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
    Q_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
    for Tx = 1:size(SIG,3)
        tmp =  hilbert(squeeze(SIG(:,:,Tx)));
        I_SIG(:,:,Tx) = real(tmp);
        Q_SIG(:,:,Tx) = imag(tmp);
    end
    [ImageTissue, ~, ~, ~, ~, ~, ~, ~, ~, ImageBone, ~, ~, ~,...
    ~, ~, ~, ~, ~, PeriParab, ~] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
        Setup, PROBE_PARAM.XR, PROBE_PARAM.ZR, PROBE_PARAM.XS, PROBE_PARAM.ZS,...
        X, Z, I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
        ReceiveAngAtPix_rad, PeriParabIn, EndoParabIn);

    PERIBOLA_COEFS = PeriParab(1:3);
    fit_curve_Periosteum = polyval(PeriParab(1:3),X);

    SegmIndPeri = PeriParab(4:end); % Raw Segmentation

    MaskPeri = ones(size(SegmIndPeri));
    MaskPeri(SegmIndPeri==-1) = NaN;
    SegmIndPeri(SegmIndPeri==-1) = 1;
    raw_peri_segmentation_Z = Z(SegmIndPeri).*MaskPeri*1e3;
    raw_peri_segmentation_X = X.*MaskPeri*1e3;

    FREQ_Transducteur = PROBE_PARAM.FREQ_Transducteur*1e6;
    wavelength_skin = C_TISSUE/(FREQ_Transducteur);
    wavelength_bone = RECON_PARAM.C_RADIAL/(FREQ_Transducteur);
%     FullImage = stitching_TBM(ImageTissue,ImageBone,ImageMarrow,P_Periosteum,P_Endosteum,X,Z,wavelength_skin,wavelength_bone);

    if plot1
        FullImage = stitching_TB(ImageTissue,ImageBone,PERIBOLA_COEFS,X,Z,wavelength_skin,wavelength_bone);
       
        close all
        FullImage = FullImage/max(FullImage(:)); 
        FullImage_dB = 20*log10(FullImage);
        fig_periosteum_segm=figure(100);
        fig_periosteum_segm.Position=[1 -153 960 1108];
        assignin('base','fig_periosteum_segm',fig_periosteum_segm);
        imagesc(X*1e3,Z*1e3,FullImage_dB)
        xticks(-20:2:20)
        colorbar
        hold on
        plot(raw_peri_segmentation_X,raw_peri_segmentation_Z,'c-','linewidth',2)
        plot(X*1e3,fit_curve_Periosteum*1e3,'r-.','linewidth',2)
        % plot(X*1e3,fit_curve_Endosteum*1e3,'r-.','linewidth',2)
        xlabel('\bf Lateral Position [mm]', 'Interpreter','latex')
        ylabel('\bf Depth [mm]', 'Interpreter','latex')
        caxis([-displayed_dynamic_range_dB 0])
        axis image
        colormap gray
        
        title({'Fine segmentation of the',' periosteal interface','Pixel size = $\frac{\lambda}{12}$'}, 'Interpreter','latex');
        set(gca,'fontsize',14)
        legend('Raw segmentation','Parabolic fit')
    end 
 
    
    % Autofocus
    HalfOpeningAngInSkinDeg=90;
    HalfOpeningAngInLensRad = asin(PROBE_PARAM.C_LENS/C_TISSUE*sind(HalfOpeningAngInSkinDeg)); % [rad]
    
    %
    % create image axes/pixels coordinates
    MIN_CORTICAL = 1.5e-3;%evalin('base','MIN_CORTICAL');
    MAX_CORTICAL = 9e-3;
    
    C_LENSBONE = C_TISSUE;
    WhichModel = 4;

    % Fill in the Setup array.
    SetupBone = [PROBE_PARAM.LENS_THICKNESS RECON_PARAM.D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL...
        HalfOpeningAngInLensRad PROBE_PARAM.Fs C_LENSBONE C_TISSUE ...
        RECON_PARAM.C_AXIAL RECON_PARAM.C_RADIAL RECON_PARAM.ANISO_SHAPE_COEF...
        RECON_PARAM.NCORE, WhichModel];
    
    X = RECON_PARAM.X; % image width [m]
    Z = RECON_PARAM.Z; % image depth [m]
  
    C_RADIAL = 3200;
    half_range = 0.10*RECON_PARAM.C_RADIAL;
    C_RADIAL_TEST_VALUES = C_RADIAL-half_range:50:C_RADIAL+half_range; % define here the range of values you want to test
    
        
        % Fill in the Setup array.
        % SetupBone = [PROBE_PARAM.LENS_THICKNESS D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL...
        %     RECON_PARAM.HalfOpeningAngInLensRad PROBE_PARAM.FS PROBE_PARAM.C_LENS ...
        %     C_TISSUE C_AXIAL C_RADIAL RECON_PARAM.ANISO_SHAPE_COEF...
        %     NCORE, WhichModel];
    tstart_autofocus_bone = datetime('now');
        [ImageTissue, ImageBone, spatial_energy_2D, ...
            sharpness_NormVariance, sharpness_Brenner, PeriParab, AutofocusResult] =...
            AutofocusBone_v4(SetupBone, PROBE_PARAM.XR, PROBE_PARAM.ZR, PROBE_PARAM.XS, PROBE_PARAM.ZS, X, Z,...
            I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, PeriParabIn, C_RADIAL_TEST_VALUES);
            
    tend_autofocus_bone = datetime('now');
        
    %
    BoneSpeed = AutofocusResult(2);
    fit_Periosteum = polyval(PeriParab,X);

    
    if plot2
        fig=figure(10);
        fig.Position= [3935 36 1755 942];
        subplot 141
        imagesc(X*1e3,Z*1e3,ImageTissue)
        hold on
        plot(X*1e3,fit_Periosteum*1e3,'r')
        xlabel('width [mm]')
        ylabel('depth [mm]')
        title('soft tissue image')
        axis image
        set(gca,'fontsize',14)
        xticks([-20:2:20])
        subplot 142
        imagesc(X*1e3,Z*1e3,ImageBone)
        xlabel('width [mm]')
        ylabel('depth [mm]')
        axis image
        set(gca,'fontsize',14)
        title('bone cortex image')
        xticks([-20:2:20])
        
        % figure(2)
        subplot 143
        plot(C_RADIAL_TEST_VALUES,spatial_energy_2D,'kv-')
        hold on
        plot(C_RADIAL_TEST_VALUES,sharpness_Brenner,'g*-')
        plot(C_RADIAL_TEST_VALUES,sharpness_NormVariance,'ms-')
        ylabel('focus quality metrics')
        xlabel('model parameter value')
        set(gca,'fontsize',14)
        legend('intensity','sharpness Brenner','sharpness NormVariance','fontsize',12, ...
            'Location','best')
        grid on
        xticks(2500:100:5000)
        
        subplot 144
        plot(C_RADIAL_TEST_VALUES,spatial_energy_2D+sharpness_Brenner+sharpness_NormVariance,'bo-')
        hold on
        plot(AutofocusResult(2),AutofocusResult(1),'rx','markersize',10,'linewidth',2)
        ylabel('sum of focus quality metrics')
        xlabel('model parameter value')
        set(gca,'fontsize',14)
        legend('sum of focus quality metrics','peak parabolic fit','fontsize',12, ...
            'Location','best')
        grid on
        xticks(2500:100:4000)
    end
end