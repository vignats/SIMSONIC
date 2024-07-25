%--- IMPORTANT ---%
% ALL PARAMETERS MUST BE OF CLASS REAL DOUBLE
% VECTORS MUST BE ROW VECTORS


%---------------------------%
%---- INPUT PARAMETERS -----%
%---------------------------%
% LENS_THICKNESS: Lens thickness (meter)
% FS: Sampling frequency (Hz)
% D_PROBE_PERIOS: a priori max depth of outer surface of bone (meter)
% MIN_CORTICAL: a priori min thickness of bone cortex (m)
% MAX_CORTICAL: a priori max thickness of bone cortex (m)

% TxHalfOpeningAngInLensRad: acceptance angle (max half opening angle) of transmit ray in lens (radian)

%---- parameter only for ImageTissue, ImageBone, ImageMarrow:
% RxHalfOpeningAngInLensRad: acceptance angle (max half opening angle) at receive ray at array element in lens (radian)


%---- parameters only for Beam_I_... and Beam_Q_... images:
% FNumberMin: fixed F-number, for large depth actual F-number might be larger
% ReceiveAngAtPix_rad: receive angles at the pixel, can be a vector for vector flow imaging [rad]
% Max_err_allowed_ReceiveAngle_rad: allowed error on receive angle at the pixel [rad]

% C_LENS: wavespeed in the lens of the probe (m/s)
% C_TISSUE: wavespeed in soft tissue between probe and bone (m/s)
% C_MARROW: wavespeed in marrow between probe and bone (m/s)
% C_AXIAL: wavespeed in bone in direction of max wavespeed (parallel to local outer surface of bone) (m/s)
% C_RADIAL: wavespeed in bone in direction of min wavespeed (m/s)
% ANISO_SHAPE_COEF: parameter of anisotropy (typically between 1 and 2), See Thomsen "Weak Elastic Anisotropy", Geophysics 1986

% ReceiveSubApertureApodis: 0, 1, 2 for rectangular, hamming, hann
% NCORE: number of CPU cores to use for OpenMP parallel computing
% ReConTo: 1, 2 or 3 (1 for soft tissue only, 2 for cutaneous tissue and cortical bone, 3 for until the marrow)
% NeedTravelTime: 0, 1 or 2 (0 for compute without saving, 1 for compute travel times/receive & transmit angles and save them, 2 for do not calculate but load travel times/receive & transmit angles)

% (XR, ZR): coordinates of array elements (meter)
% (XS, ZS): coordinates of (virtual) sources (meter)
% add_to_delay_firing: travel time from virtual sources to nearest array element (s)
% (X, Z): coordinates of image pixels (meter)
% (I_SIG, Q_SIG): in-phase and quadrature sampled RF data, 3D array dimension is (Time,Receiver index,Transmission index)
% Tx_to_be_used: indices of Sources to be used during reconstruction

% (PeriParabIn, EndoParabIn): Three parabola coefficients of periosteum and endosteum (if needed the user can force the parabolas of periosteum and/or endosteum)
% if automatic segmentation desired, must be a vector with 3 zeros:
% PeriParabIn = [ 0 0 0 ];
% EndoParabIn = [ 0 0 0 ];
% otherwise, if parabola coefficients are known, must be vector with 3 coefficients [a b c] (ax^2+bx+c)

% XMinSegmentationEndo: lateral x coordinate where the segmentation of the endosteum starts (meter)
% XMaxSegmentationEndo: lateral x coordinate where the segmentation of the endosteum ends (meter)


%---------------------------%
%---- OUTPUT PARAMETERS ----%
%---------------------------%
% ImageTissue, ImageBone, ImageMarrow: 2D image [Nz,Nx], full aperture in receive with max half opening angle
% Beam_I_... and Beam_Q_... images: 4D image [Nz,Nx,nRx,nTx], fixed receive angle at pixel and fixed F-number, for vector flow imaging
% PeriParab, EndoParab: First three elements are parabola coefficients of periosteum and endosteum,
% followed by depth indices (starting at 1 for first element of Z array, Matlab convention) of raw segmentation (-1 if no segmentation was done)


    
% Fill in the Setup array.
Setup = [LENS_THICKNESS D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL TxHalfOpeningAngInLensRad RxHalfOpeningAngInLensRad...
    FNumberMin Max_err_allowed_ReceiveAngle_rad FS C_LENS C_TISSUE C_MARROW...
    C_AXIAL C_RADIAL ANISO_SHAPE_COEF NCORE ReConTo NeedTravelTime ReceiveSubApertureApodis...
    XMinSegmentationEndo XMaxSegmentationEndo];

if ReConTo == 1
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
            Setup, XR, ZR, XS, ZS, X, Z,...
            I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
            ReceiveAngAtPix_rad, PeriParabIn, EndoParabIn);
        
elseif ReConTo == 2
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue,...
    ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone,...
    Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone,...
    APix_T_Bone, APix_R_Bone, PeriParab, EndoParab] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
        Setup, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
        ReceiveAngAtPix_rad, PeriParabIn, EndoParabIn);

elseif ReConTo == 3
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue,...
    ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone,...
    Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone, APix_T_Bone, APix_R_Bone,...
    ImageMarrow, Time_T_Marrow, Time_R_Marrow, Angle_T_Marrow,...
    Angle_R_Marrow, Beam_I_Marrow, Beam_Q_Marrow, APix_T_Marrow, APix_R_Marrow, PeriParab, EndoParab] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
        Setup, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
        ReceiveAngAtPix_rad, PeriParabIn, EndoParabIn);
    
end
        

    P_Periosteum = PeriParab(1:3);
    P_Endosteum = EndoParab(1:3);
    fit_curve_Periosteum = polyval(PeriParab(1:3),X);
    fit_curve_Endosteum = polyval(EndoParab(1:3),X);

    SegmIndPeri = PeriParab(4:end); % Raw Segmentation
    SegmIndEndo = EndoParab(4:end); % Raw Segmentation

    MaskPeri = ones(size(SegmIndPeri));
    MaskPeri(SegmIndPeri==-1) = NaN;
    SegmIndPeri(SegmIndPeri==-1) = 1;
    
    MaskEndo = ones(size(SegmIndEndo));
    MaskEndo(SegmIndEndo==-1) = NaN;
    SegmIndEndo(SegmIndEndo==-1) = 1;
    
    wavelength_skin = C_TISSUE/(FREQ_Transducteur);
    wavelength_bone = C_RADIAL/(FREQ_Transducteur);
    FullImage = stitching_TBM(ImageTissue,ImageBone,ImageMarrow,P_Periosteum,P_Endosteum,X,Z,wavelength_skin,wavelength_bone);
%     FullImage = stitching_TB(ImageTissue,ImageBone,P_Periosteum,X,Z,wavelength_skin,wavelength_bone);
    
figure(100)
imagesc(X*1e3,Z*1e3,FullImage)
hold on
plot(X.*MaskPeri*1e3,Z(SegmIndPeri).*MaskPeri*1e3,'c-','linewidth',2)
plot(X.*MaskEndo*1e3,Z(SegmIndEndo).*MaskEndo*1e3,'c-','linewidth',2)
plot(X*1e3,fit_curve_Periosteum*1e3,'r-.','linewidth',2)
plot(X*1e3,fit_curve_Endosteum*1e3,'r-.','linewidth',2)
xlabel('width [mm]')
ylabel('depth [mm]')
caxis([-displayed_dynamic_range_dB 0])
axis image
colormap gray
set(gca,'fontsize',14)
