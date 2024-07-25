clear all
close all
clc

%%----------------------------------------------------------------------------------------------%%
%%---- PROCEDURE IN THREE STEPS:
%%---- 1) RECON WITH SINGLE-ELEMENT TRANSMISSIONS FOR SEGMENTATION OF PERIOSTEUM AND ENDOSTEUM
%%---- 2) RECON FIRST FRAME WITH PLANE-WAVE TRANSMISSIONS TO CHECK BEFORE RECON ALL FRAMES
%%---- 3) RECON ALL FRAMES WITH PLANE-WAVE TRANSMISSIONS
%%----------------------------------------------------------------------------------------------%%

displayed_dynamic_range_dB = 50;
filtering_RF_data_YES = 1;
apodization_start_end_YES = 1;
TGC_YES = 1;
correction_delay_matching_layers_YES = 1;

%----------
%----------
% bone wavespeed model
C_TISSUE            = 1600; % [m/s]
C_AXIAL             = 3900; % [m/s]
C_RADIAL            = 3100; % [m/s]
ANISO_SHAPE_COEF    = 1.5; % unitless
C_MARROW            = 1420; % [m/s]
%----------
%----------

acquisition_name = 'Longi_Sujet5_Tibia_N1_HD_3.h5';

%%
%%----------------------------------------------------------------------------------------------%%
%%---- RECON WITH SINGLE-ELEMENT TRANSMISSIONS FOR SEGMENTATION OF PERIOSTEUM AND ENDOSTEUM ----%%
%%----------------------------------------------------------------------------------------------%%


file_name = ['SA_' acquisition_name];


param = h5read(file_name,'/param');
RF = h5read(file_name,'/RF');
LENGTH = param.LENGTH;
NELEMENTS = param.NELEMENTS;
CONNECTOR = param.CONNECTOR;
NSOURCES = NELEMENTS; % for single-element transmission
SIG = zeros(LENGTH,NELEMENTS,NSOURCES);

RF = RF(:,:,1); % first frame

for Tx = 1:NSOURCES
    SIG(:,:,Tx) =  RF(1 + (Tx-1)*LENGTH:Tx*LENGTH,CONNECTOR);
end
clear RF


NCORE                   = 8; % number of CPU cores for OpenMP parallel computing
ReConTo                 = 2; % 1 for soft tissue only, 2 for cutaneous tissue and cortical bone, 3 for until the marrow
NeedTravelTime          = 0; % 0, 1, 2 for compute without saving, compute travel times/receive & transmit angles and save them, do not calculate but load travel times/receive & transmit angles.

%--------------------------------------------------%
%           P4-1 probe
%--------------------------------------------------%
% lens parameters found with semblance on single reflection on interface
% silicone/water (15 oct 2019)
correction_delay    = -70.9e-9; % [s]
offset              = 9; % time to peak of the round-trip waveform [samples]
Fs                  = 10e6; % [Hz]
tone_burst_cycles   = 4;
PITCH = param.PITCH;
FREQ_Transducteur = param.FREQ_Transducteur*1e6;
C_LENS = param.C_LENS;
LENS_THICKNESS = param.LENS_THICKNESS;

%----------
%----------

% horizontal positions of the elements
XR = (0:NELEMENTS-1) *PITCH; % [m]
XR = XR - mean(XR); % place the center of coordinate system in the middle of transducer array
% vertical positions of the elements
ZR = zeros(1,length(XR)); % [m]

% coordinates of sources
XS = XR;
ZS = ZR;
NSOURCES = size(SIG,3);
% for synthetic aperture imaging, zero transmit delay added (only used for virtual sources)
add_to_delay_firing = zeros(1,NSOURCES);



% a priori knownledge on anatomy
D_PROBE_PERIOS       = 8e-3;
MIN_CORTICAL         = 2e-3;
MAX_CORTICAL         = 5e-3;

% define dimension and pixel size of reconstructed image
xmin    = -13e-3; % [m]
xmax    = 13e-3; % [m]
zmin    = 5.5e-3; % [m]
zmax    = 10.5e-3; % [m]
pixel_size = C_TISSUE/(FREQ_Transducteur*4); % [m]

% create image axes/pixels coordinates
X       = xmin:pixel_size:xmax; % image width [m]
Z       = zmin:pixel_size:zmax; % image depth [m]

% Sources to be used during reconstruction
Tx_to_be_used       = 1:NSOURCES;

% parameter only for ImageTissue, ImageBone, ImageMarrow:
HalfOpeningAngInSkinDeg = 20; % half opening angle in receive only for image reconstructed with full aperture in receive [degree]
HalfOpeningAngInLensRad = asin(C_LENS/C_TISSUE*sind(HalfOpeningAngInSkinDeg)); % [rad]
TxHalfOpeningAngInLensRad = HalfOpeningAngInLensRad;
RxHalfOpeningAngInLensRad = HalfOpeningAngInLensRad;

% Remark: no acceptance angle HalfOpeningAng used in transmit

% parameters only for Beam_I_... and Beam_Q_... images:
FnumberMin                       = 1.3; % fixed F-number, for large depth actual F-number might be larger
SubApertureApodis                = 1; % 0, 1, 2 for no windowing, hamming, hann.
% Receive angle(s)
ReceiveAngle    = deg2rad([0]); % receive angles, can be a vector for vector flow imaging [rad]
Max_err_allowed_ReceiveAngle_rad = deg2rad(1); % [rad]

% Parabolas (if needed the user can force the parabolas of periosteum and/or endosteum)
PeriParabIn = [ 0 0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros
EndoParabIn = [ 0 0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros


%%

% backward shift of raw echo signals by time to peak of the round-trip waveform
SIG = SIG(offset:end,:,:);
Nt = size(SIG,1);
Time = [0:Nt-1]/Fs;


if apodization_start_end_YES
    win = [zeros(round(2*LENS_THICKNESS/C_LENS*Fs)+offset,1)
           ones(size(SIG,1)-(round(2*LENS_THICKNESS/C_LENS*Fs)+offset),1)];
    win = repmat(win,1,size(SIG,2),size(SIG,3));
    SIG = SIG.*win;
end

if filtering_RF_data_YES
%     center_freq_bandwidth = FREQ_Transducteur;
   center_freq_bandwidth = 2.3e6;
% b = fir1(110,[0.6 1.4]*center_freq_bandwidth/(Fs/2),'bandpass'); % FullBand
b = fir1(110,[0.82 1.18]*center_freq_bandwidth/(Fs/2),'bandpass'); % NarrowBand
    for ii=1:size(SIG,2)
        for jj=1:size(SIG,3)
        SIG(:,ii,jj) = filtfilt(b,1,SIG(:,ii,jj));
        end
    end
    RF_filter = 'fir1;filtfilt;order=110;Wn=[0.82 1.18]';
else
    RF_filter = ' '; 
end

if correction_delay_matching_layers_YES
    nb_pts_fft = 2^nextpow2(Nt);
    frek = [[0:nb_pts_fft/2-1] [-nb_pts_fft/2:-1]]*Fs/nb_pts_fft;
    frek = repmat(frek.',1,NELEMENTS);
    delay_time = correction_delay*ones(1,NELEMENTS);
    for iTx=1:NSOURCES
        FFT_matrix = fft(squeeze(SIG(:,:,iTx)),nb_pts_fft);
        sig_delayed = pulse_delayingXFFT(FFT_matrix,delay_time,frek);
        SIG(:,:,iTx) = sig_delayed(1:Nt,:);
    end
end

if TGC_YES
    %%% Time Gain Compensation
    TGC = repmat((linspace(0,size(SIG,1)-1,size(SIG,1)).').^1.1,1,NELEMENTS,NSOURCES);
    TGC = TGC/max(TGC(:));
    SIG = SIG .* TGC;
end

figure(1)
subplot 131
TX_Element = 1;
Im = squeeze(SIG(:,:,TX_Element));
Im = 20*log10(abs(hilbert(Im)));
Im = Im - max(Im(:));
imagesc(XR*1e3,Time*1e6,Im)
title(['Tx with element ' int2str(TX_Element)])
xlabel('azimutal distance (mm)')
ylabel('arrival time (\mus)')
colormap gray
caxis([-displayed_dynamic_range_dB 0])
subplot 132
TX_Element = round(NSOURCES/2);
Im = squeeze(SIG(:,:,TX_Element));
Im = 20*log10(abs(hilbert(Im)));
Im = Im - max(Im(:));
imagesc(XR*1e3,Time*1e6,Im)
title(['Tx with element ' int2str(TX_Element)])
xlabel('azimutal distance (mm)')
ylabel('arrival time (\mus)')
colormap gray
caxis([-displayed_dynamic_range_dB 0])
subplot 133
TX_Element = NSOURCES;
Im = squeeze(SIG(:,:,TX_Element));
Im = 20*log10(abs(hilbert(Im)));
Im = Im - max(Im(:));
imagesc(XR*1e3,Time*1e6,Im)
title(['Tx with element ' int2str(TX_Element)])
xlabel('azimutal distance (mm)')
ylabel('arrival time (\mus)')
colormap gray
caxis([-displayed_dynamic_range_dB 0])
pause(0.5)



%%

% In-phase and quadrature sampled RF data
I_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
Q_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
for Tx = 1:size(SIG,3)
    tmp =  hilbert(squeeze(SIG(:,:,Tx)));
    I_SIG(:,:,Tx) = real(tmp);
    Q_SIG(:,:,Tx) = imag(tmp);
end



%%


XMinSegmentationEndo = X(1);
XMaxSegmentationEndo = X(end);

Setup = [LENS_THICKNESS D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL TxHalfOpeningAngInLensRad RxHalfOpeningAngInLensRad...
    FnumberMin Max_err_allowed_ReceiveAngle_rad Fs C_LENS C_TISSUE C_MARROW...
    C_AXIAL C_RADIAL ANISO_SHAPE_COEF NCORE ReConTo NeedTravelTime SubApertureApodis...
    XMinSegmentationEndo XMaxSegmentationEndo];


if ReConTo == 1
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
            Setup, XR, ZR, XS, ZS, X, Z,...
            I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
            ReceiveAngle, PeriParabIn, EndoParabIn);
        
    % sum over transmissions for image with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));

elseif ReConTo == 2
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue,...
    ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone,...
    Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone,...
    APix_T_Bone, APix_R_Bone, PeriParab, EndoParab] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
        Setup, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
        ReceiveAngle, PeriParabIn, EndoParabIn);

    P_Periosteum = PeriParab(1:3);
    P_Endosteum = EndoParab(1:3);
    fit_curve_Periosteum = polyval(PeriParab(1:3),X);
    fit_curve_Endosteum = polyval(EndoParab(1:3),X);

    % sum over transmissions for images with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));
    ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone,4).^2+sum(Beam_Q_Bone,4).^2));

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
        ReceiveAngle, PeriParabIn, EndoParabIn);
    
    P_Periosteum = PeriParab(1:3);
    P_Endosteum = EndoParab(1:3);
    fit_curve_Periosteum = polyval(PeriParab(1:3),X);
    fit_curve_Endosteum = polyval(EndoParab(1:3),X);
    
    % sum over transmissions for images with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));
    ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone,4).^2+sum(Beam_Q_Bone,4).^2));
    ImageMarrow_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Marrow,4).^2+sum(Beam_Q_Marrow,4).^2));
    ImageTissue_FixedReceiveAng = squeeze(ImageTissue_FixedReceiveAng(:,:,1));
    ImageBone_FixedReceiveAng = squeeze(ImageBone_FixedReceiveAng(:,:,1));
    ImageMarrow_FixedReceiveAng = squeeze(ImageMarrow_FixedReceiveAng(:,:,1));
end
    


%%    

%---- calculate thickness of cortical bone ----%

display_thickness_YES = 0;

if display_thickness_YES
figure(33)
plot(X*1e3,fit_curve_Periosteum*1e3,'k','linewidth',2)
hold on
plot(X*1e3,fit_curve_Endosteum*1e3,'k','linewidth',2)
axis ij
xlabel('[mm]')
ylabel('[mm]')
axis equal
end

P_midLine = mean([P_Periosteum; P_Endosteum]);
fit_curve_midLine = polyval(P_midLine,X);
if display_thickness_YES
plot(X*1e3,fit_curve_midLine*1e3,'b--','linewidth',2)
end

x1 = X;
fit_Periosteum = fit_curve_Periosteum;
for ii=[1 length(x1)]
    pt_Periosteum = [x1(ii), fit_Periosteum(ii)];
    N = [2* P_Periosteum(1)*pt_Periosteum(1) + P_Periosteum(2), -1];
    N = -N / norm(N);
    
    a = P_Endosteum(1);
    b = P_Endosteum(2);
    c = P_Endosteum(3);
    thickness1 = intersection_rayon_interface(a,b,c,pt_Periosteum,N);
    pt_Endosteum = pt_Periosteum + thickness1*N;
    if ii==1
        pt_Endosteum_left = pt_Endosteum;
    elseif ii==length(x1)
        pt_Endosteum_right = pt_Endosteum;
    end
end
X_start_Endosteum = pt_Endosteum_left(1);
X_end_Endosteum = pt_Endosteum_right(1);


x1 = X;
thicknessB = [];
thicknessT = [];
for ii=1:length(x1)
    pt_midLine = [x1(ii), fit_curve_midLine(ii)];
    N = [2* P_midLine(1)*pt_midLine(1) + P_midLine(2), -1];
    N = -N / norm(N);
    
    a = P_Endosteum(1);
    b = P_Endosteum(2);
    c = P_Endosteum(3);
    thickness = intersection_rayon_interface(a,b,c,pt_midLine,N);
    pt_Endosteum = pt_midLine + thickness*N;
    
    if (pt_Endosteum(1)>X_start_Endosteum) && (pt_Endosteum(1)<X_end_Endosteum)
        thicknessB = [thicknessB thickness];
        if display_thickness_YES
        figure(33)
        plot([pt_midLine(1) pt_Endosteum(1)]*1e3,[pt_midLine(2) pt_Endosteum(2)]*1e3,'r')
        end
        N = -N;
        a = P_Periosteum(1);
        b = P_Periosteum(2);
        c = P_Periosteum(3);
        thickness = intersection_rayon_interface(a,b,c,pt_midLine,N);
        thicknessT = [thicknessT thickness];
        pt_Periosteum = pt_midLine + thickness*N;
        if display_thickness_YES
        figure(33)
        plot([pt_midLine(1) pt_Periosteum(1)]*1e3,[pt_midLine(2) pt_Periosteum(2)]*1e3,'m')
        end
    end
end

thicknessMidLine = thicknessT+thicknessB;
mean_cortical_thickness = mean(abs(thicknessMidLine));
if display_thickness_YES
title(['mean cortical thickness = ' num2str(mean_cortical_thickness*1e3,3) ' mm']);
end
disp(['mean cortical thickness (normal to mean axis) is ' num2str(mean_cortical_thickness*1e3,3) ' mm']);

%%    


if ReConTo == 1

    figure(2)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image

elseif ReConTo == 2

    figure(2)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image

    figure(3)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image

    FullImage = ImageTissue;
    FullImage(ImageBone~=0) = ImageBone(ImageBone~=0);

    FullImage_FixedReceiveAng = ImageTissue_FixedReceiveAng;
    FullImage_FixedReceiveAng(ImageBone_FixedReceiveAng~=0) = ImageBone_FixedReceiveAng(ImageBone_FixedReceiveAng~=0);
    
    
    figure(10)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(FullImage))
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(FullImage_FixedReceiveAng))
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    
elseif ReConTo == 3

    figure(2)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image

    figure(3)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image

    figure(4)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageMarrow))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageMarrow_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    
    FullImage = ImageTissue;
    FullImage(ImageBone~=0) = ImageBone(ImageBone~=0);
    FullImage(ImageMarrow~=0) = ImageMarrow(ImageMarrow~=0);
    
    FullImage_FixedReceiveAng = ImageTissue_FixedReceiveAng;
    FullImage_FixedReceiveAng(ImageBone_FixedReceiveAng~=0) = ImageBone_FixedReceiveAng(ImageBone_FixedReceiveAng~=0);
    FullImage_FixedReceiveAng(ImageMarrow_FixedReceiveAng~=0) = ImageMarrow_FixedReceiveAng(ImageMarrow_FixedReceiveAng~=0);

    figure(10)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(FullImage))
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(FullImage_FixedReceiveAng))
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    
end

pause
close all

%%
%%--------------------------------------------------------------------------------------------%%
%%---- RECON FIRST FRAME WITH PLANE-WAVE TRANSMISSIONS TO CHECK BEFORE RECON ALL FRAMES ------%%
%%--------------------------------------------------------------------------------------------%%

file_name = ['PW_' acquisition_name];

RF=h5read(file_name,'/RF');

param=h5read(file_name,'/param');


XS = param.XS.';
ZS = param.ZS.';
add_to_delay_firing = param.add_to_delay_firing.';
LENGTH = param.LENGTH;
NELEMENTS = param.NELEMENTS;
CONNECTOR = param.CONNECTOR;
NSOURCES = param.NSOURCES;
SIG = zeros(LENGTH,NELEMENTS,NSOURCES);

RF = RF(:,:,1); % first frame

for Tx = 1:NSOURCES
    SIG(:,:,Tx) =  RF(1 + (Tx-1)*LENGTH:Tx*LENGTH,CONNECTOR);
end
clear RF



Tx_to_be_saved = 1:length(XS);
SIG = SIG(:,:,Tx_to_be_saved);
NSOURCES = length(Tx_to_be_saved);
XS = XS(Tx_to_be_saved);
ZS = ZS(Tx_to_be_saved);
add_to_delay_firing = add_to_delay_firing(Tx_to_be_saved);
% Sources to be used during reconstruction
Tx_to_be_used       = 1:NSOURCES;


% define dimension and pixel size of reconstructed image
xmin    = -8e-3; % [m]
xmax    = 8e-3; % [m]
% create image axes/pixels coordinates
X       = xmin:pixel_size:xmax; % image width [m]

ReceiveAngle    = deg2rad([-30 0 30]); % receive angles, can be a vector for vector flow imaging [rad]

NeedTravelTime                   = 1; % 0, 1, 2 for compute without saving, compute travel times/receive & transmit angles and save them, do not calculate but load travel times/receive & transmit angles.

% Parabolas (if needed the user can force the parabolas of periosteum and/or endosteum)
PeriParabIn = PeriParab(1:3); % if automatic segmentation desired, must be a vector with 5 zeros
EndoParabIn = EndoParab(1:3); % if automatic segmentation desired, must be a vector with 5 zeros

% backward shift of raw echo signals by time to peak of the round-trip waveform
SIG = SIG(offset:end,:,:);
Nt = size(SIG,1);
Time = [0:Nt-1]/Fs;

if apodization_start_end_YES
    delay_firing = zeros(NSOURCES,NELEMENTS);
    %----- transmit delays -------
    for iTx=1:NSOURCES
        delay_firing(iTx,:) = sqrt((XR-XS(iTx)).^2 + ZS(iTx).^2)/C_LENS;
        delay_firing(iTx,:) = delay_firing(iTx,:) - add_to_delay_firing(iTx);
    end
    Time_matrix = repmat(Time',1,NELEMENTS);
    muting_mask = zeros(size(SIG));
    for iTx=1:NSOURCES
        delay_firing_matrix = repmat(delay_firing(iTx,:),Nt,1);
        muting_mask(:,:,iTx) = Time_matrix > (delay_firing_matrix + 8/FREQ_Transducteur);
    end
    SIG = SIG.*muting_mask;
end

if filtering_RF_data_YES
    for ii=1:size(SIG,2)
        for jj=1:size(SIG,3)
        SIG(:,ii,jj) = filtfilt(b,1,SIG(:,ii,jj));
        end
    end
end

Nt = size(SIG,1);
Time = [0:Nt-1]/Fs;


if correction_delay_matching_layers_YES
    nb_pts_fft = 2^nextpow2(Nt);
    frek = [[0:nb_pts_fft/2-1] [-nb_pts_fft/2:-1]]*Fs/nb_pts_fft;
    frek = repmat(frek.',1,NELEMENTS);
    delay_time = correction_delay*ones(1,NELEMENTS);
    for iTx=1:NSOURCES
    FFT_matrix = fft(squeeze(SIG(:,:,iTx)),nb_pts_fft);
    sig_delayed = pulse_delayingXFFT(FFT_matrix,delay_time,frek);
    SIG(:,:,iTx) = sig_delayed(1:Nt,:);
    end
end

if TGC_YES
    %%% Time Gain Compensation
    TGC = repmat((linspace(0,size(SIG,1)-1,size(SIG,1)).').^1.1,1,NELEMENTS,NSOURCES);
    TGC = TGC/max(TGC(:));
    SIG = SIG .* TGC;
end

figure(1)
subplot 131
TX_Element = 1;
Im = squeeze(SIG(:,:,TX_Element));
Im = 20*log10(abs(hilbert(Im)));
Im = Im - max(Im(:));
imagesc(XR*1e3,Time*1e6,Im)
title(['Tx with element ' int2str(TX_Element)])
xlabel('azimutal distance (mm)')
ylabel('arrival time (\mus)')
colormap gray
caxis([-displayed_dynamic_range_dB 0])
subplot 132
TX_Element = round(NSOURCES/2);
Im = squeeze(SIG(:,:,TX_Element));
Im = 20*log10(abs(hilbert(Im)));
Im = Im - max(Im(:));
imagesc(XR*1e3,Time*1e6,Im)
title(['Tx with element ' int2str(TX_Element)])
xlabel('azimutal distance (mm)')
ylabel('arrival time (\mus)')
colormap gray
caxis([-displayed_dynamic_range_dB 0])
subplot 133
TX_Element = NSOURCES;
Im = squeeze(SIG(:,:,TX_Element));
Im = 20*log10(abs(hilbert(Im)));
Im = Im - max(Im(:));
imagesc(XR*1e3,Time*1e6,Im)
title(['Tx with element ' int2str(TX_Element)])
xlabel('azimutal distance (mm)')
ylabel('arrival time (\mus)')
colormap gray
caxis([-displayed_dynamic_range_dB 0])
pause(0.5)



%%

% In-phase and quadrature sampled RF data
I_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
Q_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
for Tx = 1:size(SIG,3)
    tmp =  hilbert(squeeze(SIG(:,:,Tx)));
    I_SIG(:,:,Tx) = real(tmp);
    Q_SIG(:,:,Tx) = imag(tmp);
end



%%


XMinSegmentationEndo = X(1);
XMaxSegmentationEndo = X(end);

Setup = [LENS_THICKNESS D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL TxHalfOpeningAngInLensRad RxHalfOpeningAngInLensRad...
    FnumberMin Max_err_allowed_ReceiveAngle_rad Fs C_LENS C_TISSUE C_MARROW...
    C_AXIAL C_RADIAL ANISO_SHAPE_COEF NCORE ReConTo NeedTravelTime SubApertureApodis...
    XMinSegmentationEndo XMaxSegmentationEndo];


if ReConTo == 1
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
            Setup, XR, ZR, XS, ZS, X, Z,...
            I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
            ReceiveAngle, PeriParabIn, EndoParabIn);
        
    % sum over transmissions for image with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));

elseif ReConTo == 2
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue,...
    ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone,...
    Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone,...
    APix_T_Bone, APix_R_Bone, PeriParab, EndoParab] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
        Setup, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
        ReceiveAngle, PeriParabIn, EndoParabIn);

    P_Periosteum = PeriParab(1:3)
    P_Endosteum = EndoParab(1:3)
    fit_curve_Periosteum = polyval(PeriParab(1:3),X);
    fit_curve_Endosteum = polyval(EndoParab(1:3),X);

    % sum over transmissions for images with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));
    ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone,4).^2+sum(Beam_Q_Bone,4).^2));
    ImageTissue_FixedReceiveAng = squeeze(ImageTissue_FixedReceiveAng(:,:,2));
    ImageBone_FixedReceiveAng = squeeze(ImageBone_FixedReceiveAng(:,:,2));
%     ImageTissue_FixedReceiveAng = squeeze(ImageTissue_FixedReceiveAng(:,:,3));
%     ImageBone_FixedReceiveAng = squeeze(ImageBone_FixedReceiveAng(:,:,3));
    
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
        ReceiveAngle, PeriParabIn, EndoParabIn);
    
    P_Periosteum = PeriParab(1:3);
    P_Endosteum = EndoParab(1:3);
    fit_curve_Periosteum = polyval(PeriParab(1:3),X);
    fit_curve_Endosteum = polyval(EndoParab(1:3),X);
    
    % sum over transmissions for images with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));
    ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone,4).^2+sum(Beam_Q_Bone,4).^2));
    ImageMarrow_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Marrow,4).^2+sum(Beam_Q_Marrow,4).^2));
    ImageTissue_FixedReceiveAng = squeeze(ImageTissue_FixedReceiveAng(:,:,2));
    ImageBone_FixedReceiveAng = squeeze(ImageBone_FixedReceiveAng(:,:,2));
    ImageMarrow_FixedReceiveAng = squeeze(ImageMarrow_FixedReceiveAng(:,:,2));
end
    
nearest_Z_Peri = min(fit_curve_Periosteum);
farthest_Z_Endo = max(fit_curve_Endosteum);
ind_Z_start = ceil(find(abs(nearest_Z_Peri-Z) == min(abs(nearest_Z_Peri-Z))));
ind_Z_end = ceil(find(abs(farthest_Z_Endo-Z) == min(abs(farthest_Z_Endo-Z))));


%%
if ReConTo == 1

    figure(2)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image

elseif ReConTo == 2

    figure(2)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image

    figure(3)
    subplot 121
    im = 20*log10(sqrt(Beam_I_Bone(:,:,1,15).^2+Beam_Q_Bone(:,:,1,15).^2));
    imagesc(X*1e3,Z*1e3,im)
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle -45 and F-number')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    im = 20*log10(sqrt(Beam_I_Bone(:,:,3,1).^2+Beam_Q_Bone(:,:,3,1).^2));
    imagesc(X*1e3,Z*1e3,im)
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle +45 and F-number')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image

    FullImage = ImageTissue;
    FullImage(ImageBone~=0) = ImageBone(ImageBone~=0);

    FullImage_FixedReceiveAng = ImageTissue_FixedReceiveAng;
    FullImage_FixedReceiveAng(ImageBone_FixedReceiveAng~=0) = ImageBone_FixedReceiveAng(ImageBone_FixedReceiveAng~=0);
    
    
    figure(10)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(FullImage))
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(FullImage_FixedReceiveAng))
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    
elseif ReConTo == 3

    figure(2)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageTissue_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    hold off
    axis image

    figure(3)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageBone_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image

    figure(4)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(ImageMarrow))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(ImageMarrow_FixedReceiveAng))
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    hold off
    axis image
    
    FullImage = ImageTissue;
    FullImage(ImageBone~=0) = ImageBone(ImageBone~=0);
    FullImage(ImageMarrow~=0) = ImageMarrow(ImageMarrow~=0);
    
    FullImage_FixedReceiveAng = ImageTissue_FixedReceiveAng;
    FullImage_FixedReceiveAng(ImageBone_FixedReceiveAng~=0) = ImageBone_FixedReceiveAng(ImageBone_FixedReceiveAng~=0);
    FullImage_FixedReceiveAng(ImageMarrow_FixedReceiveAng~=0) = ImageMarrow_FixedReceiveAng(ImageMarrow_FixedReceiveAng~=0);

    figure(10)
    subplot 121
    imagesc(X*1e3,Z*1e3,20*log10(FullImage))
    title('full aperture in receive with max half opening angle')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    subplot 122
    imagesc(X*1e3,Z*1e3,20*log10(FullImage_FixedReceiveAng))
    title('fixed receive angle and F-number')
    hold on
    plot(X*1e3,fit_curve_Periosteum*1e3,'k:')
    plot(X*1e3,fit_curve_Endosteum*1e3,'k:')
    xlabel('width [mm]')
    ylabel('depth [mm]')
    axis image
    
end


pause

%%
%%----------------------------------------------------------%%
%%---- RECON ALL FRAMES WITH PLANE-WAVE TRANSMISSIONS ------%%
%%----------------------------------------------------------%%

RF=h5read(file_name,'/RF');

FRAMERATE = param.FRAMERATE;

nb_frames = size(RF,3);

NeedTravelTime                   = 2; % 0, 1, 2 for compute without saving, compute travel times/receive & transmit angles and save them, do not calculate but load travel times/receive & transmit angles.

% Parabolas (if needed the user can force the parabolas of periosteum and/or endosteum)
PeriParabIn = [ 0 0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros
EndoParabIn = [ 0 0 0 ]; % if automatic segmentation desired, must be a vector with 5 zeros

% Tx_to_be_saved = 1:length(Tx_to_be_used);%1:2:length(XS);
Tx_Compound = zeros(3,15);
Tx_Compound(1,:) = 1:15;
Tx_Compound(2,1:5) = 1:5;
Tx_Compound(3,1:5) = 11:15;


% Memory allocation
ImageBone_FixedReceiveAng_I = zeros(length(Z),length(X),length(ReceiveAngle),size(Tx_Compound,1),nb_frames);
ImageBone_FixedReceiveAng_Q = ImageBone_FixedReceiveAng_I;

for iT=1:nb_frames
    
tmp = RF(:,:,iT);
SIG = zeros(LENGTH,NELEMENTS,NSOURCES);
for Tx = 1:NSOURCES
    SIG(:,:,Tx) =  tmp(1 + (Tx-1)*LENGTH:Tx*LENGTH,CONNECTOR);
end
clear tmp

%%

% backward shift of raw echo signals by time to peak of the round-trip waveform
SIG = SIG(offset:end,:,:);
Nt = size(SIG,1); % number of temporal samples
Time = [0:Nt-1]/Fs;


if apodization_start_end_YES
    delay_firing = zeros(NSOURCES,NELEMENTS);
    %----- transmit delays -------
    for iTx=1:NSOURCES
        delay_firing(iTx,:) = sqrt((XR-XS(iTx)).^2 + ZS(iTx).^2)/C_LENS;
        delay_firing(iTx,:) = delay_firing(iTx,:) - add_to_delay_firing(iTx);
    end
    Time_matrix = repmat(Time',1,NELEMENTS);
    muting_mask = zeros(size(SIG));
    for iTx=1:NSOURCES
        delay_firing_matrix = repmat(delay_firing(iTx,:),Nt,1);
        muting_mask(:,:,iTx) = Time_matrix > (delay_firing_matrix + 8/FREQ_Transducteur);
    end
    SIG = SIG.*muting_mask;
end


if filtering_RF_data_YES
    for ii=1:size(SIG,2)
        for jj=1:size(SIG,3)
        SIG(:,ii,jj) = filtfilt(b,1,SIG(:,ii,jj));
        end
    end
end

Nt = size(SIG,1);
Time = [0:Nt-1]/Fs;


if correction_delay_matching_layers_YES
    nb_pts_fft = 2^nextpow2(Nt);
    frek = [[0:nb_pts_fft/2-1] [-nb_pts_fft/2:-1]]*Fs/nb_pts_fft;
    frek = repmat(frek.',1,NELEMENTS);
    delay_time = correction_delay*ones(1,NELEMENTS);
    for iTx=1:NSOURCES
    FFT_matrix = fft(squeeze(SIG(:,:,iTx)),nb_pts_fft);
    sig_delayed = pulse_delayingXFFT(FFT_matrix,delay_time,frek);
    SIG(:,:,iTx) = sig_delayed(1:Nt,:);
    end
end

if TGC_YES
    %%% Time Gain Compensation
    TGC = repmat((linspace(0,size(SIG,1)-1,size(SIG,1)).').^1.1,1,NELEMENTS,NSOURCES);
    TGC = TGC/max(TGC(:));
    SIG = SIG .* TGC;
end


%%

% In-phase and quadrature sampled RF data
I_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
Q_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
for Tx = 1:size(SIG,3)
    tmp =  hilbert(squeeze(SIG(:,:,Tx)));
    I_SIG(:,:,Tx) = real(tmp);
    Q_SIG(:,:,Tx) = imag(tmp);
end



%%


XMinSegmentationEndo = X(1);
XMaxSegmentationEndo = X(end);

Setup = [LENS_THICKNESS D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL TxHalfOpeningAngInLensRad RxHalfOpeningAngInLensRad...
    FnumberMin Max_err_allowed_ReceiveAngle_rad Fs C_LENS C_TISSUE C_MARROW...
    C_AXIAL C_RADIAL ANISO_SHAPE_COEF NCORE ReConTo NeedTravelTime SubApertureApodis...
    XMinSegmentationEndo XMaxSegmentationEndo];


if ReConTo == 1
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
            Setup, XR, ZR, XS, ZS, X, Z,...
            I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
            ReceiveAngle, PeriParabIn, EndoParabIn);
        
    % sum over transmissions for image with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue,4).^2+sum(Beam_Q_Tissue,4).^2));

elseif ReConTo == 2
    
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue,...
    ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone,...
    Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone,...
    APix_T_Bone, APix_R_Bone, PeriParab, EndoParab] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
        Setup, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
        ReceiveAngle, PeriParabIn, EndoParabIn);

    P_Periosteum = PeriParab(1:3);
    P_Endosteum = EndoParab(1:3);
    fit_curve_Periosteum = polyval(PeriParab(1:3),X);
    fit_curve_Endosteum = polyval(EndoParab(1:3),X);

    % sum over transmissions for images with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageBone_FixedReceiveAng_I(:,:,:,1,iT) = sum(Beam_I_Bone(:,:,:,Tx_Compound(1,:)),4); % transmit compounding
    ImageBone_FixedReceiveAng_Q(:,:,:,1,iT) = sum(Beam_Q_Bone(:,:,:,Tx_Compound(1,:)),4); % transmit compounding
    ImageBone_FixedReceiveAng_I(:,:,:,2,iT) = sum(Beam_I_Bone(:,:,:,Tx_Compound(2,1:5)),4); % transmit compounding
    ImageBone_FixedReceiveAng_Q(:,:,:,2,iT) = sum(Beam_Q_Bone(:,:,:,Tx_Compound(2,1:5)),4); % transmit compounding
    ImageBone_FixedReceiveAng_I(:,:,:,3,iT) = sum(Beam_I_Bone(:,:,:,Tx_Compound(3,1:5)),4); % transmit compounding
    ImageBone_FixedReceiveAng_Q(:,:,:,3,iT) = sum(Beam_Q_Bone(:,:,:,Tx_Compound(3,1:5)),4); % transmit compounding
    
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
        ReceiveAngle, PeriParabIn, EndoParabIn);
    
    P_Periosteum = PeriParab(1:3);
    P_Endosteum = EndoParab(1:3);
    fit_curve_Periosteum = polyval(PeriParab(1:3),X);
    fit_curve_Endosteum = polyval(EndoParab(1:3),X);
    
    % sum over transmissions for images with fixed receive angle and F-number
    % if multiple receive angles, then receive angle must be chosen with
    % dimension 3 in Beam_I and Beam_Q
    ImageBone_FixedReceiveAng_I(:,:,:,:,iT) = sum(Beam_I_Bone,4); % transmit compounding
    ImageBone_FixedReceiveAng_Q(:,:,:,:,iT) = sum(Beam_Q_Bone,4); % transmit compounding
    
end
    
iT

end


ImageBone_FixedReceiveAng_I = ImageBone_FixedReceiveAng_I(ind_Z_start:ind_Z_end,:,:,:,:);
ImageBone_FixedReceiveAng_Q = ImageBone_FixedReceiveAng_Q(ind_Z_start:ind_Z_end,:,:,:,:);
Z = Z(ind_Z_start:ind_Z_end);
APix_T_Bone = APix_T_Bone(:,ind_Z_start:ind_Z_end,:);
APix_R_Bone = APix_R_Bone(:,ind_Z_start:ind_Z_end,:);

save([file_name(1:end-3) '_recon_for_LSVectorFlow_Bone.mat'],'ImageBone_FixedReceiveAng_I',...
    'ImageBone_FixedReceiveAng_Q','P_Periosteum','P_Endosteum','Tx_to_be_used','D_PROBE_PERIOS','MIN_CORTICAL','MAX_CORTICAL','RF_filter',...
    'FnumberMin','SubApertureApodis','fit_curve_Periosteum','fit_curve_Endosteum','FRAMERATE','X','Z','ReceiveAngle','TGC_YES',...
    'APix_T_Bone','APix_R_Bone','Tx_Compound','C_TISSUE','C_AXIAL','C_RADIAL','ANISO_SHAPE_COEF','C_MARROW','center_freq_bandwidth','-v7.3')


