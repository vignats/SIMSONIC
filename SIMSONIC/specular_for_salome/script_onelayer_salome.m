% ----------------------
% Add directories that contain relevant functions
addpath(genpath("~/CODE"))
addpath(genpath('~/Softwares/'))
% -----------------------
close all force
tstart = datetime("now");
clc
fprintf('----------------------------------------------------\n')
fprintf('Started the processing at %s\n', tstart)
% Directory that contains simulation data
simul_data_dir = '~/1D_SPECULAR_MODEL';%'/calculSSD/adia/SIMULATION/SPECULAR/FLUID';
%
cortical_porosity = evalin('base','cortical_porosity');
pore_size = evalin('base','pore_size');
transmit_sequence = evalin('base','transmit_sequence'); 
geometry = lower(evalin('base','geometry'));
clearvars -except pixel_size_wl xmin xmax zmin zmax ...
    simul_data_dir cortical_porosity pore_size transmit_sequence geometry tstart
data_dir = sprintf('%s/PORE_SIZE_%02d',simul_data_dir,pore_size);
path_to_rf = sprintf('%s/%s_recorded_rf_no_atten_%s_porosity_%02d_pore_size_%d.mat', ...
    data_dir,transmit_sequence, geometry,cortical_porosity, pore_size);

recorded  = load(path_to_rf);

%
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

[Angle_T,Angle_R] = function_get_angle_of_view_sa(X,Z,PROBE_PARAM.XS,PROBE_PARAM.ZS,...
    PROBE_PARAM.XR,PROBE_PARAM.ZR);
[Time_T, Time_R] = function_get_time_of_flight_sa(X,Z,BEST_C_LAYER_1,PROBE_PARAM.XS,PROBE_PARAM.ZS, ...
    PROBE_PARAM.XR,PROBE_PARAM.ZR);

TIME_of_FLIGHT.Time_T =Time_T;
TIME_of_FLIGHT.Time_R = Time_R;
ANGLE_of_VIEW.Angle_T = Angle_T;
ANGLE_of_VIEW.Angle_R = Angle_R;
%% Delayed signal for each pixel with combination of each transmit/receive
fprintf('----------Full map of delayed signals ----------\n');
tic
FULL_DELAYED_MAP_OF_RF_DATA_PIXEL = function_get_angular_rf_interpol(RECON_PARAM,TIME_of_FLIGHT, ...
        PROBE_PARAM.FS);
toc

%% IMAGE RECONSTRUCTION with DAS algorithm
DAS_IMAGE = zeros([numel(Z) numel(X)]);
% tx_fnumber = .7;
% rx_fnumber = .7;
half_receive_ap_tiss = atand(1/2/rx_fnumber);
half_receive_ap_bone = atand(1/2/rx_fnumber);
half_transmi_ap = atand(1/2/tx_fnumber);
for iz=1:numel(Z)
    for ix=1:numel(X)
        delayed_signals_tiss = FULL_DELAYED_MAP_OF_RF_DATA_PIXEL(:,:,iz,ix);
        delayed_signals_tiss(isnan(delayed_signals_tiss))=0;
        at_tiss_px = ANGLE_of_VIEW.Angle_T(:,iz,ix);
        ar_tiss_px = ANGLE_of_VIEW.Angle_R(:,iz,ix);
        t_idx_tiss = find(at_tiss_px>=-half_transmi_ap & at_tiss_px<=half_transmi_ap);
        r_idx_tiss = find(ar_tiss_px>=-half_receive_ap_tiss & ar_tiss_px<=half_receive_ap_tiss);
        low_rez_tiss = sum(tukeywin(numel(r_idx_tiss),1).'.*delayed_signals_tiss(:,r_idx_tiss),2);
        DAS_IMAGE(iz,ix) = abs(sum((tukeywin(numel(t_idx_tiss),0).*low_rez_tiss(t_idx_tiss))));%abs(sum((delayed_signals_tiss(t_idx_tiss,r_idx_tiss))));%
    end
end
%

DAS_IMAGE(isnan(DAS_IMAGE))=0;
DAS_IMAGE = DAS_IMAGE/max(DAS_IMAGE(:));


DAS_IMAGE_dB = 20*log10(DAS_IMAGE/max(DAS_IMAGE(:)));
% ------------------------
% Affichage de l'image reconstruite en échelle logarithmique
figure, pcolor(Xmm, Zmm,DAS_IMAGE_dB), axis ij image, shading flat, 
xlabel('Lateral position [mm]', Interpreter='latex')
ylabel('Depth [mm]', Interpreter='latex')
caxis([-40 0])
colormap gray

%% SPECULAR TRANSFORM
tilt_max = 70; NTilts = 2*tilt_max+1;
TILT_angles = linspace(-tilt_max,tilt_max,NTilts);

fprintf('----------Getting 1D specular transform----------\n');
fprintf('\tSoft Tissues\n')
tic
[SPECULAR_TRANSFORM] = function_get_specular_transform_interpol(RECON_PARAM,...
    TIME_of_FLIGHT,ANGLE_of_VIEW, ANGLE_of_VIEW, PROBE_PARAM.FS, TILT_angles);
toc
%% Illustration
ix = find(X>0e-3,1);
iz = find(Z>=5e-3,1);

s=-1;
ar_tiss_px = ANGLE_of_VIEW.Angle_R(:,iz,ix);
at_tiss_px = ANGLE_of_VIEW.Angle_T(:,iz,ix);
ar_tiss_elem = ANGLE_of_VIEW.Angle_R(:,iz,ix);
at_tiss_elem = ANGLE_of_VIEW.Angle_T(:,iz,ix);
directivity_map_tiss = ones([96 96]);
% directivity_map_tiss(:,abs(at_tiss_elem)>atand(1/2/tx_fnumber))=0;
% directivity_map_tiss(abs(ar_tiss_elem)>atand(1/2/rx_fnumber),:)=0;

delayed_signals_perios = s*(FULL_DELAYED_MAP_OF_RF_DATA_PIXEL(:,:,iz,ix)).*directivity_map_tiss;
cplx_spe_transform_perios = SPECULAR_TRANSFORM(:,iz,ix);
spe_transform_perios = s*real(SPECULAR_TRANSFORM(:,iz,ix));
das_transform = sum(s*real(delayed_signals_perios),1);

fig_illust=figure(Position=[1317 1 1244 1321]) ;
tilt = 0;
tlo=tiledlayout(2,1);
nexttile(tlo,1,[1 1])
pcolor(at_tiss_px,ar_tiss_px, real(delayed_signals_perios))
colorbar
hold on 
plot(at_tiss_px,-at_tiss_px+2*(tilt),'--')
xticks(-100:20:100)
yticks(-100:20:100)
shading flat 
axis image
axis([at_tiss_px(end) at_tiss_px(1) at_tiss_px(end) at_tiss_px(1)])
xlabel('\bf Transmit angle $\gamma_t [^{\circ}]$','Interpreter','latex')
ylabel('\bf Receive angle $\gamma_r [^{\circ}]$','Interpreter','latex')
title({'Periosteum','Delayed signals '},interpreter='latex')

nexttile(tlo,2)
hold on
plot(TILT_angles,spe_transform_perios)%/max(spe_transform_perios))
xline(tilt,'--')
xticks(-100:20:100)
axis tight
xlabel('\bf Specular orientation $\beta [^{\circ}]$','Interpreter','latex')
% ylabel('\bf Depth [mm]','Interpreter','latex')
title('Specular Transform',interpreter='latex')
%% MODEL OF SPECULAR TRANSFORM : simplified version
model_version = 'exact_1';%exact_0 exact_1 simplified_0 simplified_1
excitation_signal=hilbert(SimSonic2DReadSgl([simul_data_dir '/Signal.sgl']));

PROBE_PARAM.XR = PROBE_PARAM.XS;
fprintf('---------GET SPECULAR TRANSFORM MODELS---------\n')

fprintf('\tSoft Tissues\n')
tic
% MODEL_MAP_TISSUE = get_simplified_model_1D_tissu_NEW(estimated_geometry,PROBE_PARAM, ...
%     TIME_of_FLIGHT,excitation_signal,PROBE_PARAM.FS,ANGLE_of_VIEW,TILT_angles);
NZ = numel(Z); NX = numel(X);
MODEL_MAP = zeros([NTilts NZ NX]);
SPECULAR_INTERFACE = [0 0]; % coefficient d'une droite [a b] z=ax+b
wb=waitbar(0,'Simplified specuar model...');
for iz=1:NZ
    waitbar(iz/NZ,wb,'Simplified specuar model ...');
    zp = Z(iz);
    for ix=1:NX
        transmit_time = TIME_of_FLIGHT.Time_T(:,iz,ix);
        transmit_angle = ANGLE_of_VIEW.Angle_R(:,iz,ix);
        xp = X(ix);
        if all(transmit_time)
            SPECULAR_MODEL=function_get_pix_simplified_tissue_model(xp,zp,...
                SPECULAR_INTERFACE,BEST_C_LAYER_1,PROBE_PARAM.XS, ...
                excitation_signal,PROBE_PARAM.FS);
            % Interpolation to match desired specular angles
            SPECULAR_MODEL = interp1(transmit_angle,SPECULAR_MODEL,TILT_angles.');
            MODEL_MAP(:,iz,ix) = SPECULAR_MODEL;
        end
    end
end
toc
close(wb)
%% MAP OF SPECULARITY
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

thresh = .5;
D_THETA = (TILT_angles(2)-TILT_angles(1));
lags = D_THETA*(-NTilts+1:NTilts-1);
[PROBA_MAP_TISS, ind_T] = max(abs(CORRELATION),[],3);
PROBA_MAP_TISS(isnan(PROBA_MAP_TISS))=0;


SPECULAR_TILT_MAP_TISS = lags(ind_T);
PROBA_MAP_TISS_bin  = double((PROBA_MAP_TISS>thresh));
PROBA_MAP_TISS_bin(PROBA_MAP_TISS_bin==0)=nan;
TILT_MAP_TISS_bin = SPECULAR_TILT_MAP_TISS.*PROBA_MAP_TISS_bin;

%
ST = SPECULAR_TILT_MAP_TISS;
ST(isnan(SPECULAR_TILT_MAP_TISS)) = 0;

figure('Position', [1 1 1920 1108]),
tlo = tiledlayout(1,2);
nexttile(tlo,1)
pcolor(Xmm,Zmm,PROBA_MAP_TISS)
xlabel('Lateral position [mm]', Interpreter='latex')
ylabel('Depth [mm]', Interpreter='latex')
colorbar
title('Probability map TISSUE',interpreter='latex')
axis ij image, shading flat
clim([0 1])


nexttile(tlo,2)
pcolor(Xmm,Zmm,SPECULAR_TILT_MAP_TISS)
colorbar
title('Detected specular angle TISSUE',interpreter='latex')
% colormap gray
clim([-1 1]*45)
axis ij image, shading flat
xlabel('Lateral position [mm]', Interpreter='latex')
ylabel('Depth [mm]', Interpreter='latex')