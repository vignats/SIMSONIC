addpath /home_local/Experimental/BoneImageVersionGuillaume/PostProcessing_v4/'Autofocus SoftTissue_Jan2023'/
addpath /home_local/Experimental/BoneImageVersionGuillaume/PostProcessing_v4/'Autofocus Bone_Jan2023'/
C_TISSUE_TEST_VALUES = C_TISSUE-150:10:C_TISSUE+150; % define here the range of values you want to test    
xmaxsoft = 7e-3; zmaxsoft = 8e-3;
PITCH = 300e-6;%33*9e-6;%
XR = (0:NRx-1)*PITCH;XR = XR-mean(XR);ZR = zeros(size(XR));
XSOFT = -xmaxsoft:pixel_size:xmaxsoft;
ZSOFT = zmin:pixel_size:zmaxsoft;
LENS_THICKNESS      = 1e-9;%1.3e-3; % [m]

HalfOpeningAngInSkinDeg=90;
HalfOpeningAngInLensRad = asin(C_LENS/C_TISSUE*sind(HalfOpeningAngInSkinDeg)); % [rad]

Setup_soft = [LENS_THICKNESS ...
HalfOpeningAngInLensRad PROBE_PARAM.Fs C_LENS C_TISSUE NCORE];

Tx_to_be_used = 1:NTx;
add_to_delay_firing = zeros([1 NTx]);
tic
[ImageTissue, spatial_energy_2D, sharpness_NormVariance, sharpness_Brenner, AutofocusResult] = ...
    AutofocusSoftTissue_v4(Setup_soft, XR, ZR, XR, ZR, XSOFT, ZSOFT,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, C_TISSUE_TEST_VALUES);
toc

%
% Show metrics
figure('Position', [1921 102 1280 708]),
subplot 131
pcolor(XSOFT*1e3,ZSOFT*1e3,ImageTissue)
ax= gca;
ax.YDir = 'reverse';
shading interp
xlabel('width [mm]')
ylabel('depth [mm]')
title('soft tissue image')
axis image
set(gca,'fontsize',14)
xticks(-20:2:20)

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

BEST_C_TISSSUE = AutofocusResult(2);


f=gcf;
print(f,[save_dir '/../../AUTOFOCUS_' simul_type '_' sim_name...
'.png'],'-dpng');...,'-fillpage')
print(f,[save_dir '/../../AUTOFOCUS_' simul_type '_' sim_name...
'.pdf'],'-dpdf','-fillpage')
%%
xmaxbone = 10e-3; zmaxbone= 12e-3;
fc = 1;
XBONE = -xmaxbone:pixel_size/fc:xmaxbone;
zminbone = 2e-3;
ZBONE= zminbone:pixel_size/fc:zmaxbone;

C_RADIAL = 3100;
PeriParabIn = [0 0 0];%[30.7533    15.4919e-3   5.387e-3];%
WhichModel = 4; 
C_RADIAL_TEST_VALUES = C_RADIAL-600:50:C_RADIAL+600; % define here the range of values you want to test

MIN_CORTICAL = 2e-3;
MAX_CORTICAL = 6e-3;
C_LENSBONE = BEST_C_TISSSUE;
% Fill in the Setup array.
SetupBone = [LENS_THICKNESS D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL...
    HalfOpeningAngInLensRad PROBE_PARAM.FS C_LENSBONE BEST_C_TISSSUE ...
    C_AXIAL C_RADIAL ANISO_SHAPE_COEF...
    NCORE, WhichModel];
tic
    [ImageTissue, ImageBone, spatial_energy_2D, sharpness_NormVariance, ...
        sharpness_Brenner, PeriParab, AutofocusResult] =...
        AutofocusBone_v4(SetupBone, XR, ZR, XR, ZR, XBONE, ZBONE,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, PeriParabIn, C_RADIAL_TEST_VALUES);
        
toc
BEST_C_RADIAL = AutofocusResult(2);

%%
fit_Periosteum = polyval(PeriParab,XBONE);

fig=figure(10);
fig.Position= [3935 36 1755 942];
subplot 141
imagesc(XBONE*1e3,ZBONE*1e3,ImageTissue)
hold on
plot(XBONE*1e3,fit_Periosteum*1e3,'r')
xlabel('width [mm]')
ylabel('depth [mm]')
title('soft tissue image')
axis image
set(gca,'fontsize',14)
xticks(-20:2:20)
subplot 142
pcolor(XBONE*1e3,ZBONE*1e3,ImageBone)
ax = gca;ax.YDir = "reverse";
shading flat
xlabel('width [mm]')
ylabel('depth [mm]')
axis image
set(gca,'fontsize',14)
title('bone cortex image')
xticks(-20:2:20)

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
f=gcf;
print(f,[save_dir '/../../AUTOFOCUSBONE_' simul_type '_' sim_name...
'.png'],'-dpng');...,'-fillpage')
print(f,[save_dir '/../../AUTOFOCUSBONE_' simul_type '_' sim_name...
'.pdf'],'-dpdf','-fillpage')


%%
C_TISSUE = BEST_C_TISSSUE;C_RADIAL = BEST_C_RADIAL;