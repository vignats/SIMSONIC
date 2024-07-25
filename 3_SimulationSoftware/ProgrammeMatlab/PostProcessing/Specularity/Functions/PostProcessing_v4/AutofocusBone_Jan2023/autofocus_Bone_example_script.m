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
% D_PROBE_PERIOS: a priori max depth of outer surface of bone (meter)
% MIN_CORTICAL: a priori min thickness of bone
% MAX_CORTICAL: a priori max thickness of bone
% HalfOpeningAngInLensRad: max half opening angle at receive element (radian)
% C_TISSUE: known wavespeed in soft tissue between probe and bone (m/s)
% C_AXIAL: known (or guess of) wavespeed in bone in direction of max wavespeed (m/s)
% C_RADIAL: known (or guess of) wavespeed in bone in direction of min wavespeed (m/s)
% ANISO_SHAPE_COEF: known (or guess of) anisotropy shape factor (unitless)
% NCORE: number of CPU cores to use for OMP parallel computing
% WhichModel: 1, 2, 3 or 4 (see below)
% model_Values: row vector containing parameter values you want to test
% (X_El, Z_El): coordinates of array elements (meter)
% (XS, ZS): coordinates of (virtual) sources (meter)
% add_to_delay_firing: travel time from virtual sources to nearest array element (second), zeros for single-element transmissions
% (X, Z): coordinates of image pixels (meter)
% (I_SIG, Q_SIG): in-phase and quadrature sampled RF data, 3D array dimension is (Time,Receiver index,Transmission index)
% Tx_to_be_used: indices of transmissions to be used

% (PeriParabIn): Three parabola coefficients of periosteum (if needed the user can force the parabolas of periosteum)
% if automatic segmentation desired, must be a vector with 3 zeros:
% PeriParabIn = [ 0 0 0 ];
% otherwise, if parabola coefficients are known, must be vector with 3 coefficients [a b c] (ax^2+bx+c)


    
%---- Select which parameter to test ----%
%     1 = CAxial (Anisotropic velocity model, CRadial and AnisoShape coef are known)
%     2 = CRadial (Anisotropic velocity model, CAxial and AnisoShape coef are known)
%     3 = AnisoShape coef (Anisotropic velocity model, CAxial and CRadial are known)
%     4 = CRadial (Isotropic velocity model)
%----------------------------------------%
    WhichModel = 4; 
    model_Values = C_RADIAL-250:50:C_RADIAL+250; % define here the range of values you want to test

    
    % Fill in the Setup array.
    Setup = [LENS_THICKNESS D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL...
        HalfOpeningAngInLensRad FS C_LENS C_TISSUE C_AXIAL C_RADIAL ANISO_SHAPE_COEF...
        NCORE, WhichModel];

    [ImageTissue, ImageBone, spatial_energy_2D, sharpness_NormVariance, sharpness_Brenner, PeriParab, AutofocusResult] =...
        AutofocusBone_v4(Setup, X_El, Z_El, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, PeriParabIn, model_Values);
        

    
%%

fit_Periosteum = polyval(PeriParab,X);

figure(1)
subplot 121
imagesc(X*1e3,Z*1e3,ImageTissue)
hold on
plot(X*1e3,fit_Periosteum*1e3,'r')
xlabel('width [mm]')
ylabel('depth [mm]')
title('soft tissue image')
axis image
set(gca,'fontsize',14)
subplot 122
imagesc(X*1e3,Z*1e3,ImageBone)
xlabel('width [mm]')
ylabel('depth [mm]')
axis image
set(gca,'fontsize',14)
title('bone cortex image')

figure(2)
subplot 121
plot(model_Values,spatial_energy_2D,'kv-')
hold on
plot(model_Values,sharpness_Brenner,'g*-')
plot(model_Values,sharpness_NormVariance,'ms-')
ylabel('focus quality metrics')
xlabel('model parameter value')
set(gca,'fontsize',14)
legend('intensity','sharpness Brenner','sharpness NormVariance','fontsize',12)
grid on
subplot 122
plot(model_Values,spatial_energy_2D+sharpness_Brenner+sharpness_NormVariance,'bo-')
hold on
plot(AutofocusResult(2),AutofocusResult(1),'rx','markersize',10,'linewidth',2)
ylabel('sum of focus quality metrics')
xlabel('model parameter value')
set(gca,'fontsize',14)
legend('sum of focus quality metrics','peak parabolic fit','fontsize',12)
grid on

