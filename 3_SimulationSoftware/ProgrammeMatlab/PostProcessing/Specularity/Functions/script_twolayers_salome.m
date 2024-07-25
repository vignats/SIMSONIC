% ----------------------
% Add directories that contain relevant functions
% addpath(genpath("~/CODE"))
% addpath(genpath('~/Softwares/'))
% -----------------------
close all force
tstart = datetime("now");
clc
fprintf('----------------------------------------------------\n')
fprintf('Started the processing at %s\n', tstart)
% Directory that contains simulation data
simul_data_dir = '~/1D_SPECULAR_MODEL';%'/calculSSD/adia/SIMULATION/SPECULAR/FLUID';
%

% data_dir = sprintf('%s/PORE_SIZE_%02d',simul_data_dir,pore_size);
data_dir =  '/calculSSD/salome/Simulation-Test/simulation_two'
parameters = load(fullfile(data_dir, 'parameters.mat'));
recorded  = LoadRfData(parameters.probe, data_dir);

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
zmin = 2e-3;%
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

%%
normalize_bone_image = 1;
XMIN_ENDO = xmin; XMAX_ENDO=xmax;
tic
[DAS_output,PERI,ENDO, Angles, Times] = DAS_boneimage_Salome(RECON_PARAM,rx_fnumber, ...
XMIN_ENDO,XMAX_ENDO);
toc

DAS_image_dB = 20*log10(DAS_output);
% close all
figure,
imagesc(Xmm,Zmm, DAS_image_dB)
axis ij image
xlabel('Lateral position [mm]')
ylabel('Depth [mm]')

% TIME OF FLIGHT and ANGLES OF VIEW
RECON_PARAM.PERI = PERI;
clear RAW_RF_SIG RAW_IQ_SIG RAW_Q_SIG RAW_I_SIG 
XT=PROBE_PARAM.XS;ZT = PROBE_PARAM.ZS; % SA transmission

TOF_TISS.Time_T = Times.Time_T_Tissue;
TOF_TISS.Time_R = Times.Time_R_Tissue;
AOV_px_deg_TISS.Angle_T = rad2deg(Angles.APix_T_Tissue);
AOV_px_deg_TISS.Angle_R = rad2deg(Angles.APix_T_Tissue);
TOF_BONE.Time_T = Times.Time_T_Bone;
TOF_BONE.Time_R = Times.Time_R_Bone;
AOV_px_deg_BONE.Angle_T = rad2deg(Angles.APix_T_Bone);
AOV_px_deg_BONE.Angle_R = rad2deg(Angles.APix_T_Bone);
AOV_elem_deg_BONE.Angle_T = rad2deg(Angles.Angle_T_Bone);
AOV_elem_deg_BONE.Angle_R = rad2deg(Angles.Angle_T_Bone);
AOV_elem_deg_TISS.Angle_T = rad2deg(Angles.Angle_T_Tissue);
AOV_elem_deg_TISS.Angle_R = rad2deg(Angles.Angle_R_Tissue);
%% Get RF map wrt angles 
fprintf('----------Getting RF MAP----------\n');
fprintf('\tSoft Tissues\n')
    tic
    RF_MAP_ANGULAR_TISS = function_get_angular_rf_interpol(RECON_PARAM,TOF_TISS, PROBE_PARAM.FS);
    toc

fprintf('\tBone Tissues\n')
    tic
    RF_MAP_ANGULAR_BONE = function_get_angular_rf_interpol(RECON_PARAM,TOF_BONE, PROBE_PARAM.FS);
    toc

%% Get specular transform
tilt_max = 70; NTilts = 2*tilt_max+1;
TILT_angles = linspace(-tilt_max,tilt_max,NTilts);

fprintf('----------Getting 1D specular transform----------\n');
fprintf('\tSoft Tissues\n')
tic
[SPECULAR_TRANSFORM_TISS] = function_get_specular_transform_interpol(RECON_PARAM,...
    TOF_TISS,AOV_px_deg_TISS, AOV_elem_deg_TISS, PROBE_PARAM.FS, TILT_angles);
toc
%
fprintf('\tBone Tissues\n')
tic
[SPECULAR_TRANSFORM_BONE] = function_get_specular_transform_interpol(RECON_PARAM,...
    TOF_BONE,AOV_px_deg_BONE,AOV_elem_deg_BONE, PROBE_PARAM.FS, TILT_angles);
toc


%% SPECULAR TRANSFORM MODELS
model_version = 'exact_1';%exact_0 exact_1 simplified_0 simplified_1
excitation_signal=hilbert(SimSonic2DReadSgl([data_dir '/Signal.sgl']));
PROBE_PARAM.XR = PROBE_PARAM.XS;
fprintf('---------GET SPECULAR TRANSFORM MODELS---------\n')
estimated_geometry = struct('X',X,'Z',Z,'lens_thick',PROBE_PARAM.LENS_THICKNESS,...
 ...
    'fit_curve_Periosteum',Z(PERI(4:end)),...
    'peribola_coef',PERI(1:3),'endobola_coef',ENDO(1:3),...
    'V0',C_LENS, ...
    'V1',RECON_PARAM.C_TISSUE, ...
    'V2',RECON_PARAM.C_RADIAL...
    );
fprintf('\tSoft Tissues\n')
tic
MODEL_MAP_TISSUE = get_simplified_model_1D_tissu_NEW(estimated_geometry,PROBE_PARAM, ...
    TOF_TISS,excitation_signal,PROBE_PARAM.FS,AOV_px_deg_TISS,TILT_angles);
% get_model_map_v1_tissue_simul(model_version, ...
%     estimated_geometry,PROBE_PARAM.XR, TILT_angles,...
%     AOV_px_deg_TISS,AOV_elem_deg_TISS,TOF_TISS,excitation_signal,PROBE_PARAM.FS);
toc
MODEL_MAP_TISSUE(isnan(MODEL_MAP_TISSUE)) = 0;

fprintf('\tBone Tissues - 1D model\n')
tic,
MODEL_MAP_BONE= get_simplified_model_1D_NEW(estimated_geometry,PROBE_PARAM, ...
    TOF_BONE,excitation_signal,PROBE_PARAM.FS,AOV_px_deg_BONE,TILT_angles);
toc
MODEL_MAP_BONE(isnan(MODEL_MAP_BONE)) = 0;


%% MAP OF SPECULARITY
tic
CORR3_BONE_nested= xcorr_nested(SPECULAR_TRANSFORM_BONE, ...
    MODEL_MAP_BONE);
toc
tic
CORR3_TISS_nested = xcorr_nested(SPECULAR_TRANSFORM_TISS, ...
    MODEL_MAP_TISSUE);
toc
%
%

thresh = .5;
D_THETA = (TILT_angles(2)-TILT_angles(1));
lags = D_THETA*(-NTilts+1:NTilts-1);
[PROBA_MAP_TISS, ind_T] = max(abs(CORR3_TISS_nested),[],3);
PROBA_MAP_TISS(isnan(PROBA_MAP_TISS))=0;
[PROBA_MAP_BONE, ind_B] = max(abs(CORR3_BONE_nested),[],3);
PROBA_MAP_BONE(isnan(PROBA_MAP_BONE))=0;

SPECULAR_TILT_MAP_TISS = lags(ind_T);
PROBA_MAP_TISS_bin  = double((PROBA_MAP_TISS>thresh));
PROBA_MAP_TISS_bin(PROBA_MAP_TISS_bin==0)=nan;
TILT_MAP_TISS_bin = SPECULAR_TILT_MAP_TISS.*PROBA_MAP_TISS_bin;

SPECULAR_TILT_MAP_BONE = lags(ind_B);
PROBA_MAP_BONE_bin  = double((PROBA_MAP_BONE>thresh));
PROBA_MAP_BONE_bin(PROBA_MAP_BONE_bin==0)=nan;
TILT_MAP_BONE_bin = SPECULAR_TILT_MAP_BONE.*PROBA_MAP_BONE_bin;

%
ST = SPECULAR_TILT_MAP_TISS;
ST(isnan(SPECULAR_TILT_MAP_TISS)) = 0;

figure('Position', [1 1 1920 1108]),
tlo = tiledlayout(1,2);
nexttile(tlo,1)
pcolor(Xmm,Zmm,PROBA_MAP_BONE)
xlabel('Lateral position [mm]', 'Interpreter', 'latex')
ylabel('Depth [mm]', 'Interpreter', 'latex')
colorbar
title('Probability map Bone', 'Interpreter','latex')
axis ij image, shading flat
clim([0 1])

nexttile(tlo,2)
pcolor(Xmm,Zmm,SPECULAR_TILT_MAP_BONE)
colorbar
title('Detected specular angle Bone','Interpreter','latex')
% colormap gray
clim([-1 1]*45)
axis ij image, shading flat
xlabel('Lateral position [mm]', 'Interpreter','latex')
ylabel('Depth [mm]', 'Interpreter','latex')
%%