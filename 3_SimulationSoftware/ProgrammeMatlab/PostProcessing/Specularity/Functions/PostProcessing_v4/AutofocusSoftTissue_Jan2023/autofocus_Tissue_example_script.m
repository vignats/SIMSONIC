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
% (I_SIG, Q_SIG): in-phase and quadrature sampled RF data, 3D array dimension is (Time,Receiver index,Transmission index)
% Tx_to_be_used: indices of transmissions to be used
% C_TISSUE_TEST_VALUES: row vector containing values you want to test


C_TISSUE_TEST_VALUES = c_tissue-150:25:c_tissue+150; % define here the range of values you want to test

% Fill in the Setup array.
Setup = [LENS_THICKNESS HalfOpeningAngInLensRad FS C_LENS C_TISSUE NCORE];

[ImageTissue, spatial_energy_2D, sharpness_NormVariance, sharpness_Brenner, AutofocusResult] = AutofocusSoftTissue_v4(Setup, X_El, Z_El, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing, C_TISSUE_TEST_VALUES);

%%

figure(1)
imagesc(X*1e3,Z*1e3,ImageTissue)
xlabel('width [mm]')
ylabel('depth [mm]')
title('soft tissue image')
axis image
set(gca,'fontsize',14)

figure(2)
subplot 121
plot(C_TISSUE_TEST_VALUES,spatial_energy_2D,'kv-')
hold on
plot(C_TISSUE_TEST_VALUES,sharpness_Brenner,'g*-')
plot(C_TISSUE_TEST_VALUES,sharpness_NormVariance,'ms-')
ylabel('focus quality metrics')
xlabel('model parameter value')
set(gca,'fontsize',14)
legend('intensity','sharpness Brenner','sharpness NormVariance','fontsize',12)
grid on
subplot 122
plot(C_TISSUE_TEST_VALUES,spatial_energy_2D+sharpness_Brenner+sharpness_NormVariance,'bo-')
hold on
plot(AutofocusResult(2),AutofocusResult(1),'rx','markersize',10,'linewidth',2)
ylabel('sum of focus quality metrics')
xlabel('model parameter value')
set(gca,'fontsize',14)
legend('sum of focus quality metrics','peak parabolic fit','fontsize',12)
grid on

