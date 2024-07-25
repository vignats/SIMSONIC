% ----------------------
% Add directories that contain relevant functions
% addpath(genpath("~/CODE"))
% addpath(genpath('~/Softwares/'))
% -----------------------
close all force
clearvars
tstart = datetime("now");
clc
fprintf('----------------------------------------------------\n')
fprintf('Started the processing at %s\n', tstart)
% Directory that contains simulation data
simul_data_dir = '~/1D_SPECULAR_MODEL';%'/calculSSD/adia/SIMULATION/SPECULAR/FLUID';


path_to_rf = "/calculSSD/salome/Simulation-10mai/Bone227G-Image2002"

recorded  = load(path_to_rf);

%%
C_LENS = 1540; LENS_THICKNESS = 1e-6;

RAW_RF_SIG = recorded.Signals;
PROBE_PARAM.Nsamples = recorded.NbTimePoints; % Number of time samples
PROBE_PARAM.Fs = 1/recorded.Temporal_step_us*1e6; % in Hz
PROBE_PARAM.FS = PROBE_PARAM.Fs;
PROBE_PARAM.PITCH = recorded.Pitch*recorded.Spatial_step_mm*1e-3;
PROBE_PARAM.NELEMENTS = recorded.NbElements; % Number of receivers
PROBE_PARAM.NSOURCES = size(RAW_RF_SIG,3);   % Number of emissions
PROBE_PARAM.FREQ_Transducteur = 2.5e6;
PROBE_PARAM.ELEM_WIDTH = recorded.Width*recorded.Spatial_step_mm*1e-3;

PROBE_PARAM.XS = (0:PROBE_PARAM.NSOURCES-1)*PROBE_PARAM.PITCH;
PROBE_PARAM.XS = PROBE_PARAM.XS - mean(PROBE_PARAM.XS);


%
pixel_size_wl = 1/4; % 1/4 de longueur d'onde

xmin = -10e-3;%
xmax = -xmin;%
zmin = 0.5e-3;%
zmax = 15e-3;%

BEST_C_LAYER_1 = 1540;
BEST_C_LAYER_2 = 3500;

pixel_size = pixel_size_wl*BEST_C_LAYER_1/PROBE_PARAM.FREQ_Transducteur;
X = xmin:pixel_size:xmax;
Z = zmin:pixel_size:zmax;
Xmm = X*1e3; Zmm = Z*1e3;



%

PROBE_PARAM.ZS = zeros(size(PROBE_PARAM.XS));
PROBE_PARAM.C_LENS = C_LENS;
PROBE_PARAM.LENS_THICKNESS = LENS_THICKNESS;
field_probe_struct = {'FS','FREQ_Transducteur','PITCH',...
'C_LENS','LENS_THICKNESS','NELEMENTS','Nsamples',...
'XS','ZS','NSOURCES','offset'};
nb_wavelengths = 3;
input_signal_duration = nb_wavelengths/PROBE_PARAM.FREQ_Transducteur;
PROBE_PARAM.offset = round(input_signal_duration/2*PROBE_PARAM.Fs);
assert(all(isfield(PROBE_PARAM,field_probe_struct)))

%



%
field_recon_struc = {
    'X','Z','MIN_CORTICAL', 'MAX_CORTICAL','D_PROBE_PERIOS',...
    'C_TISSUE','C_RADIAL','C_AXIAL','C_MARROW','ANISO_SHAPE_COEF',...
    'FnumberMin', 'SubApertureApodis', 'NeedTravelTime','ReConTo',...
    'Max_err_allowed_ReceiveAngle_rad','HalfOpeningAngInLensRad'};
D_PROBE_PERIOS = 10e-3;
tx_fnumber = 0.6; % Transmit fnumber
rx_fnumber = 0.6; %  Receive fnumber
RECON_PARAM = struct('MIN_CORTICAL',3e-3, 'MAX_CORTICAL',15e-3,...
    'D_PROBE_PERIOS',D_PROBE_PERIOS,...
    'C_TISSUE',BEST_C_LAYER_1,'C_RADIAL',BEST_C_LAYER_2,'C_AXIAL',BEST_C_LAYER_2,...
    'C_MARROW',1540,'ANISO_SHAPE_COEF',0,...
    'FnumberMin',rx_fnumber, 'SubApertureApodis',2, ...
    'NeedTravelTime',0,'ReConTo',2,...
    'Max_err_allowed_ReceiveAngle_rad', deg2rad(5));
TxHalfOpeningAngInSkinDeg=atand(1/2/tx_fnumber);
RECON_PARAM.HalfOpeningAngInLensRad = asin(PROBE_PARAM.C_LENS/RECON_PARAM.C_TISSUE*sind(TxHalfOpeningAngInSkinDeg)); % [rad]



RECON_PARAM.X = X; RECON_PARAM.Z = Z;
assert(all(isfield(RECON_PARAM,field_recon_struc)))
NCORE = 24;
RECON_PARAM.rf_data = RAW_RF_SIG;
RECON_PARAM.PROBE_PARAM = PROBE_PARAM;
RECON_PARAM.Setup = [RECON_PARAM.ANISO_SHAPE_COEF LENS_THICKNESS PROBE_PARAM.PITCH...
    RECON_PARAM.D_PROBE_PERIOS RECON_PARAM.MIN_CORTICAL RECON_PARAM.MAX_CORTICAL ...
    RECON_PARAM.HalfOpeningAngInLensRad rx_fnumber...
    RECON_PARAM.Max_err_allowed_ReceiveAngle_rad PROBE_PARAM.FS ...
    PROBE_PARAM.C_LENS RECON_PARAM.C_TISSUE RECON_PARAM.C_MARROW...
    RECON_PARAM.C_AXIAL RECON_PARAM.C_RADIAL NCORE ...
    RECON_PARAM.ReConTo RECON_PARAM.NeedTravelTime RECON_PARAM.SubApertureApodis];


PROBE_PARAM.XR = PROBE_PARAM.XS;
PROBE_PARAM.ZR = PROBE_PARAM.ZS;
SIG = RECON_PARAM.rf_data(PROBE_PARAM.offset+1:end,:,:);
NSOURCES = size(SIG,3);
    Tx_to_be_used       = 1:NSOURCES;

add_to_delay_firing = zeros(1,NSOURCES);
% I/Q separation.
I_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
Q_SIG = ones(size(SIG,1),size(SIG,2),size(SIG,3));
for Tx = 1:size(SIG,3)
    tmp =  hilbert(squeeze(SIG(:,:,Tx)));
    I_SIG(:,:,Tx) = real(tmp);
    Q_SIG(:,:,Tx) = imag(tmp);
end

%% Autofocus soft tissues
%
skip=1;
if skip
    BEST_C_LAYER_1 = 1540;
    RECON_PARAM.C_TISSUE = BEST_C_LAYER_1;
else
    %--- IMPORTANT ---%
    % ALL PARAMETERS MUST BE OF CLASS REAL DOUBLE
    % VECTORS MUST BE ROW VECTORS
    
    
    %--- WARNING ---%
    % the autofocus approach requires a rather larger number of different
    % transmit beams (different insonofication angles), typically at least 20
    
    
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
    half_range = 0.10*RECON_PARAM.C_TISSUE;
    step = 0.01*RECON_PARAM.C_TISSUE;
    C_TISSUE = RECON_PARAM.C_TISSUE;
    C_TISSUE_TEST_VALUES = C_TISSUE-half_range:25:C_TISSUE+half_range; % define here the range of values you want to test    
    
    % Fill in the Setup array.
    % Setup = [LENS_THICKNESS HalfOpeningAngInLensRad FS C_LENS C_TISSUE NCORE];
    Setup_soft = [PROBE_PARAM.LENS_THICKNESS ...
        RECON_PARAM.HalfOpeningAngInLensRad PROBE_PARAM.Fs PROBE_PARAM.C_LENS C_TISSUE NCORE];
    tstart_autofocus_soft = datetime('now');
    
    [ImageTissue, spatial_energy_2D, sharpness_NormVariance, sharpness_Brenner, AutofocusResult] = ...
        AutofocusSoftTissue_v4(Setup_soft, PROBE_PARAM.XR, PROBE_PARAM.ZR, ...
        PROBE_PARAM.XS, PROBE_PARAM.ZS, RECON_PARAM.X, RECON_PARAM.Z,...
            I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, C_TISSUE_TEST_VALUES);
    tend_autofocus_soft = datetime('now');
    BEST_C_LAYER_1 = AutofocusResult(2);
    RECON_PARAM.C_TISSUE = BEST_C_LAYER_1;
    combined_metric = spatial_energy_2D+sharpness_Brenner+sharpness_NormVariance;
    % Show metrics
    figure('Position', [1921 102 1280 708]),
    subplot 131
    imagesc(X*1e3,Z*1e3,ImageTissue)
    xlabel('width [mm]')
    ylabel('depth [mm]')
    title('soft tissue image')
    axis image
    set(gca,'fontsize',14)
    xticks([-20:2:20])
    
    % figure(2)
    subplot 132
    plot(C_TISSUE_TEST_VALUES,spatial_energy_2D,'kv-')
    hold on
    plot(C_TISSUE_TEST_VALUES,sharpness_Brenner,'g*-')
    plot(C_TISSUE_TEST_VALUES,sharpness_NormVariance,'ms-')
    ylabel('focus quality metrics')
    xlabel('model parameter value')
    set(gca,'fontsize',14)
    legend('intensity','sharpness Brenner','sharpness NormVariance', ...
        'fontsize',12, 'Location', 'best')
    grid on
    xticks(1000:50:2000)
    
    subplot 133
    plot(C_TISSUE_TEST_VALUES,spatial_energy_2D+sharpness_Brenner+sharpness_NormVariance,'bo-')
    hold on
    plot(AutofocusResult(2),AutofocusResult(1),'rx','markersize',10,'linewidth',2)
    ylabel('sum of focus quality metrics')
    xlabel('model parameter value')
    set(gca,'fontsize',14)
    legend('sum of focus quality metrics','peak parabolic fit','fontsize',12 ...
        ,'Location', 'best')
    grid on
    xticks(1000:50:2000)
end
%% Fine recon for periosteum segmentation

C_TISSUE = RECON_PARAM.C_TISSUE;
 

ReceiveAngAtPix_rad = deg2rad([0]); % receive angles, can be a vector for vector flow imaging [rad]
RxHalfOpeningAngInLensRad = RECON_PARAM.HalfOpeningAngInLensRad;
TxHalfOpeningAngInLensRad = deg2rad(90);
fine_resolution = C_TISSUE/(PROBE_PARAM.FREQ_Transducteur*12); % [m]
FS = PROBE_PARAM.FS;
C_LENS = PROBE_PARAM.C_LENS;
XR = PROBE_PARAM.XR;ZR = PROBE_PARAM.ZR;
XS =XR;ZS=ZR;
XMinSegmentationEndo = xmin;
XMaxSegmentationEndo = xmax;
%
MIN_CORTICAL = RECON_PARAM.MIN_CORTICAL;%
% create image axes/pixels coordinates
z1 = 17e-3;MAX_CORTICAL=3e-3;
zmin = evalin('base','zmin');
zmax    = evalin('base','zmax');%8e-3;%max(z1,D_PROBE_PERIOS+MAX_CORTICAL); % [m]
X       = xmin:fine_resolution:xmax; % image width [m]
Z       = zmin:fine_resolution:zmax; % image depth [m]

ReConTo = 2;  
C_AXIAL = C_TISSUE;
C_RADIAL = C_TISSUE;
C_MARROW = C_TISSUE;
LENS_THICKNESS = PROBE_PARAM.LENS_THICKNESS;
% Fill in the Setup array.
% Setup = [LENS_THICKNESS D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL TxHalfOpeningAngInLensRad RxHalfOpeningAngInLensRad...
%     FNumberMin Max_err_allowed_ReceiveAngle_rad FS C_LENS C_TISSUE C_MARROW...
%     C_AXIAL C_RADIAL ANISO_SHAPE_COEF NCORE ReConTo NeedTravelTime ReceiveSubApertureApodis...
%     XMinSegmentationEndo XMaxSegmentationEndo];

Setup = [PROBE_PARAM.LENS_THICKNESS RECON_PARAM.D_PROBE_PERIOS...
    RECON_PARAM.MIN_CORTICAL RECON_PARAM.MAX_CORTICAL ...
    RECON_PARAM.HalfOpeningAngInLensRad RECON_PARAM.HalfOpeningAngInLensRad...
    RECON_PARAM.FnumberMin RECON_PARAM.Max_err_allowed_ReceiveAngle_rad PROBE_PARAM.FS ...
    PROBE_PARAM.C_LENS C_TISSUE RECON_PARAM.C_MARROW...
    C_AXIAL C_RADIAL RECON_PARAM.ANISO_SHAPE_COEF NCORE...
    ReConTo RECON_PARAM.NeedTravelTime RECON_PARAM.SubApertureApodis...
    XMinSegmentationEndo XMaxSegmentationEndo];
PeriParabIn = [0 0 0];EndoParabIn=[0 0 0];
tstart_fine_recon = datetime('now');
    [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue,...
    Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue,...
    ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone,...
    Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone,...
    APix_T_Bone, APix_R_Bone, PeriParab, EndoParab] =...
    PostProcessingFNumber_WeakAniso_SUPERFAST_RxApodAll_v4(...
        Setup, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing,...
        ReceiveAngAtPix_rad, PeriParabIn, EndoParabIn);
tend_fine_recon = datetime('now');
%   
displayed_dynamic_range_dB = 40;

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
    wavelength_bone = C_RADIAL/(FREQ_Transducteur);
%     FullImage = stitching_TBM(ImageTissue,ImageBone,ImageMarrow,P_Periosteum,P_Endosteum,X,Z,wavelength_skin,wavelength_bone);
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


%% Autofocus bone tissues

    
%---- Select which parameter to test ----%
%     1 = CAxial (Anisotropic velocity model, CRadial and AnisoShape coef are known)
%     2 = CRadial (Anisotropic velocity model, CAxial and AnisoShape coef are known)
%     3 = AnisoShape coef (Anisotropic velocity model, CAxial and CRadial are known)
%     4 = CRadial (Isotropic velocity model)
%----------------------------------------%
C_TISSUE = RECON_PARAM.C_TISSUE;
PeriParabIn = PERIBOLA_COEFS;
HalfOpeningAngInSkinDeg=90;
HalfOpeningAngInLensRad = asin(PROBE_PARAM.C_LENS/C_TISSUE*sind(HalfOpeningAngInSkinDeg)); % [rad]

%
% create image axes/pixels coordinates
MIN_CORTICAL = 1.5e-3;%evalin('base','MIN_CORTICAL');
MAX_CORTICAL=9e-3;
WhichModel = 4; 
C_LENSBONE = C_TISSUE;
% Fill in the Setup array.
SetupBone = [LENS_THICKNESS D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL...
    HalfOpeningAngInLensRad PROBE_PARAM.FS C_LENSBONE C_TISSUE ...
    RECON_PARAM.C_AXIAL RECON_PARAM.C_RADIAL RECON_PARAM.ANISO_SHAPE_COEF...
    NCORE, WhichModel];

X       = xmin:pixel_size:xmax; % image width [m]
Z       = zmin:pixel_size:zmax; % image depth [m]



    
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
        AutofocusBone_v4(SetupBone, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, PeriParabIn, C_RADIAL_TEST_VALUES);
        
tend_autofocus_bone = datetime('now');
    
%

fit_Periosteum = polyval(PeriParab,X);

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

%% Elapsed time
tend = datetime('now');
fprintf('--------------------####--------------------\n')
fprintf('All Elapsed time: %s\n',tend-tstart)
fprintf('Autofocus soft tissues : %s\n', tend_autofocus_soft-tstart_autofocus_soft)
fprintf('Fine recon perisoteum  : %s\n', tend_fine_recon-tstart_fine_recon)
fprintf('Autofocus bone tissues : %s\n', tend_autofocus_bone-tstart_autofocus_bone)
fprintf('--------------------####--------------------\n')


