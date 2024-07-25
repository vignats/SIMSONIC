clear all
close all
clc

mean_filter_YES = 0
SVD_filter_YES = 0
Trad_HP_filter_YES = 1

load Longi_Sujet5_Tibia_N1_HD_3_recon_for_LSVectorFlow_Bone

nb_Tx = size(ImageBone_FixedReceiveAng_I,4);
% nb_Rx = length(ReceiveAngle);

ReceiveAngle_to_be_used = 1:3;
nb_Rx = length(ReceiveAngle_to_be_used);
ImageBone_FixedReceiveAng_I = ImageBone_FixedReceiveAng_I(:,:,ReceiveAngle_to_be_used,:,:);
ImageBone_FixedReceiveAng_Q = ImageBone_FixedReceiveAng_Q(:,:,ReceiveAngle_to_be_used,:,:);

index_zero_RxAng = find(ReceiveAngle==0);

APix_T_Bone_to_be_used = zeros(size(Tx_Compound,1),length(Z),length(X));
APix_T_Bone_to_be_used(1,:,:) = mean(APix_T_Bone(Tx_Compound(1,:),:,:),1); % transmit compounding !!!!!!!!!
APix_T_Bone_to_be_used(2,:,:) = mean(APix_T_Bone(Tx_Compound(2,1:5),:,:),1); % transmit compounding !!!!!!!!!
APix_T_Bone_to_be_used(3,:,:) = mean(APix_T_Bone(Tx_Compound(3,1:5),:,:),1); % transmit compounding !!!!!!!!!

APix_R_Bone_to_be_used = APix_R_Bone(ReceiveAngle_to_be_used,:,:);



%%
%----------------------%
%--- Clutter Filter ---%
%----------------------%
nb_frames = size(ImageBone_FixedReceiveAng_I,5);

if mean_filter_YES
    ImageBone_FixedReceiveAng_I = ImageBone_FixedReceiveAng_I - repmat(mean(ImageBone_FixedReceiveAng_I,5),1,1,1,1,size(ImageBone_FixedReceiveAng_I,5));
    ImageBone_FixedReceiveAng_Q = ImageBone_FixedReceiveAng_Q - repmat(mean(ImageBone_FixedReceiveAng_Q,5),1,1,1,1,size(ImageBone_FixedReceiveAng_Q,5));
end

if SVD_filter_YES

tre=[5 50]; % SVD
tmpI = zeros(size(ImageBone_FixedReceiveAng_I));
tmpQ = zeros(size(ImageBone_FixedReceiveAng_Q));
%---- Bone ----%
for iTx=1:nb_Tx
disp(['SVD clutter filter Bone - Tx ' int2str(iTx) '/' int2str(nb_Tx)])    
    for iRx=1:nb_Rx
        IQ_im(:,:,1,1,:)=ImageBone_FixedReceiveAng_I(:,:,iRx,iTx,:) + 1j*ImageBone_FixedReceiveAng_Q(:,:,iRx,iTx,:);
        cas=get_casorati_from_acquisition(real(IQ_im),imag(IQ_im),nb_frames,1,size(ImageBone_FixedReceiveAng_I,2),:,1,1);
        tmp=filterWithSVD(cas.ImI,size(IQ_im,5),size(IQ_im,2),size(IQ_im,1),tre);
        tmpI(:,:,iRx,iTx,:) = tmp.frames;
        tmp=filterWithSVD(cas.ImQ,size(IQ_im,5),size(IQ_im,2),size(IQ_im,1),tre);
        tmpQ(:,:,iRx,iTx,:) = tmp.frames;
    end
end
ImageBone_FixedReceiveAng_I = tmpI;
ImageBone_FixedReceiveAng_Q = tmpQ;
clear tmpI tmpQ
% clutter_filter_type = 'SVD-Threshold=10';
% function_to_apply_clutter_filter = 'filterWithSVD';

end


if Trad_HP_filter_YES

N = 5;
Wn = [0.05];
[B,A] = butter(N,Wn,'high');
% Wn = [0.05 0.25];
% [B,A] = butter(N,Wn,'bandpass');

tmpI = zeros(size(ImageBone_FixedReceiveAng_I));
tmpQ = zeros(size(ImageBone_FixedReceiveAng_Q));
%---- Bone ----%
for ix=1:length(X)
disp(['clutter filter Bone - X pixel ' int2str(ix) '/' int2str(length(X))])
    for iz=1:length(Z)
        for iRx=1:nb_Rx
            for iTx=1:nb_Tx
                tmp = squeeze(ImageBone_FixedReceiveAng_I(iz,ix,iRx,iTx,:));
                tmp = filtfilt(B,A,tmp);
%                 tmp = filter(B,A,tmp);
                tmpI(iz,ix,iRx,iTx,:) = tmp;
                tmp = squeeze(ImageBone_FixedReceiveAng_Q(iz,ix,iRx,iTx,:));
                tmp = filtfilt(B,A,tmp);
%                 tmp = filter(B,A,tmp);
                tmpQ(iz,ix,iRx,iTx,:) = tmp;
            end
        end
    end
end
ImageBone_FixedReceiveAng_I = tmpI;
ImageBone_FixedReceiveAng_Q = tmpQ;
clear tmpI tmpQ
% clutter_filter_type = 'Butterworth-N=5-Wn=0.08';
% function_to_apply_clutter_filter = 'filtfilt';

end


disp('clutter filter finished')


%%
%-----------------------------%
%------ POWER DOPPLER --------%
%-----------------------------%

kernel_size_spatial_smoothing_x = 1;
kernel_size_spatial_smoothing_z = 1;
Temporal_smooting_kernel_size = 40; % 1 does nothing, 0 average all frames, >1 temporal smoothing

PowerDoppler = zeros(size(ImageBone_FixedReceiveAng_I));
PowerDoppler = permute(PowerDoppler,[3 4 1 2 5]);
PowerDoppler_SpatialAv = zeros(size(ImageBone_FixedReceiveAng_I,3),size(ImageBone_FixedReceiveAng_I,4),size(ImageBone_FixedReceiveAng_I,5));
for iRx=1:nb_Rx
    for iTx=1:nb_Tx
        tmp = squeeze(ImageBone_FixedReceiveAng_I(:,:,iRx,iTx,:)+1j*ImageBone_FixedReceiveAng_Q(:,:,iRx,iTx,:));
        [PowerDoppler(iRx,iTx,:,:,:),PowerDoppler_SpatialAv(iRx,iTx,:)] = compute_PowerDoppler(tmp,kernel_size_spatial_smoothing_x,kernel_size_spatial_smoothing_z,Temporal_smooting_kernel_size);
    end
end
PowerDoppler = permute(PowerDoppler,[3 4 1 2 5]);


PowerDoppler_SpatialAv_ALLTxRx0 = squeeze(mean(PowerDoppler_SpatialAv,[1,2]));

start_frame = 50;
figure(33)
plot(start_frame:nb_frames,PowerDoppler_SpatialAv_ALLTxRx0(start_frame:end),'r')
ylabel({'power doppler','spatial average',['temporal smoothing (' int2str(Temporal_smooting_kernel_size) ')']})
xlabel('frame number')


ZeroRxAngInd = find(ReceiveAngle==0);

PowerDoppler_ALLTxRx = squeeze(mean(PowerDoppler(:,:,ZeroRxAngInd,:,:),4)); % !! 0 deg receive ang

disp('power doppler complete')

%%
%-----------------------------%
%--------- LS - BONE ---------%
%-----------------------------%

Spatial_smooting_kernel_size_x = 5;
Spatial_smooting_kernel_size_z = 5;

% selected_iTx = [1 1 1 2 2 2 3 3 3];
% selected_iRx = [1 2 3 1 2 3 1 2 3];
selected_iTx = [1 1];
selected_iRx = [1 3];
TxAng1 = APix_T_Bone_to_be_used(selected_iTx(1),:,:);
TxAng1 = mean(TxAng1(:),'omitnan')
TxAng2 = APix_T_Bone_to_be_used(selected_iTx(2),:,:);
TxAng2 = mean(TxAng2(:),'omitnan')
RxAng1 = APix_R_Bone_to_be_used(selected_iRx(1),:,:);
RxAng1 = mean(RxAng1(:),'omitnan')
RxAng2 = APix_R_Bone_to_be_used(selected_iRx(2),:,:);
RxAng2 = mean(RxAng2(:),'omitnan')
PSF_angle1_deg = rad2deg(mean([TxAng1 RxAng1]))
PSF_angle2_deg = rad2deg(mean([TxAng2 RxAng2]))

[vz_LS_Bone0,vx_LS_Bone0,LS_sd_residuals,LS_sd_error_Vz0,LS_sd_error_Vx0]=compute_VFlow_temporal_smooth_Bone_residual_selectTxRxpairs(ImageBone_FixedReceiveAng_I,ImageBone_FixedReceiveAng_Q,APix_T_Bone_to_be_used,APix_R_Bone_to_be_used,C_RADIAL,C_AXIAL,ANISO_SHAPE_COEF,center_freq_bandwidth,FRAMERATE,Spatial_smooting_kernel_size_x,Spatial_smooting_kernel_size_z,Temporal_smooting_kernel_size,selected_iTx,selected_iRx);


disp('LS vector flow finished')


%%

%-----------------------------%
%--------- TO - BONE ---------%
%-----------------------------%
ImIQ = squeeze(mean(ImageBone_FixedReceiveAng_I(:,:,ZeroRxAngInd,:,:) + 1j*ImageBone_FixedReceiveAng_Q(:,:,ZeroRxAngInd,:,:),4));

param.fps=FRAMERATE;
param.Nx=Spatial_smooting_kernel_size_x;  % spatial smoothing
param.Nz=Spatial_smooting_kernel_size_z;  % spatial smoothing
param.Nt=Temporal_smooting_kernel_size; % temporal smoothing
param.lag=1;
param.display=1;
param.X = X; % [m]
param.Z = Z; % [m]
param.mask_ROI = ones(length(Z),length(X));
param.assumed_spatial_period_image = C_RADIAL/center_freq_bandwidth/2; % [m]

% TO parameters
param.sigmax = 1.5; % [mm] inverse width gaussian mask in Fourier domain
param.lambdax = 1.5;  % [mm] inverse of distance between two Gaussians in Fourier domain (kx,kz)
TO_PSF_angle_deg = asind(C_RADIAL/center_freq_bandwidth/2*1000/param.lambdax)

[Vz_Full_Bone0,vz_TO_Bone0,vx_TO_Bone0,Im1,Im2]=TO_Dec2022(ImIQ,param);

disp('TO vector flow finished')


%%
%---------------------------%
%--------- DISPLAY ---------%
%---------------------------%


NT=10;
ht=ones(1,NT)/NT;
vx_LS_Bone=FiltFiltM(ht,1,vx_LS_Bone0,3);
vz_LS_Bone=FiltFiltM(ht,1,vz_LS_Bone0,3);
vx_TO_Bone=FiltFiltM(ht,1,vx_TO_Bone0,3)/1000;
vz_TO_Bone=FiltFiltM(ht,1,Vz_Full_Bone0,3)/1000;
PowerDoppler_SpatialAv_ALLTxRx = FiltFiltM(ht,1,PowerDoppler_SpatialAv_ALLTxRx0,1);



VELO=sqrt(vz_LS_Bone.^2+vx_LS_Bone.^2);
ANGL=atan2(vx_LS_Bone,vz_LS_Bone);

VELO_TO=sqrt(vz_TO_Bone.^2+vx_TO_Bone.^2);
ANGL_TO=atan2(vx_TO_Bone,vz_TO_Bone);


PD_dB_Threshold = -40;

PowerDoppler_dB = 20*log10(PowerDoppler_ALLTxRx);
tmp = PowerDoppler_dB(:,:,150:end);
PowerDoppler_dB = PowerDoppler_dB - max(tmp(:));
PD_threshold_mask = (PowerDoppler_dB>PD_dB_Threshold);
PD_threshold_mask = double(PD_threshold_mask);
PD_threshold_mask(PD_threshold_mask==0) = NaN;
PD_threshold_mask = PD_threshold_mask(:,:,1:end-1);


%%


rect = [50 50 450 580]; %[left bottom width height];
H1=figure(111);
set(111,'Position',rect)

under_samp = 2;
scaling_factor=1/30;


index_range_display = [40:2:2000];
PD = PowerDoppler_SpatialAv_ALLTxRx(index_range_display(1):index_range_display(end))/max(PowerDoppler_SpatialAv_ALLTxRx(index_range_display(1):index_range_display(end)));
time = [0:length(PD)-1]/FRAMERATE;
t_ind = 1;
clear F
for time_index = index_range_display
figure(111)
subplot(5,1,[1:2])
s2=imagesc(X*1e3,Z*1e3,PowerDoppler_dB(:,:,time_index));
caxis([PD_dB_Threshold 0])
hold on
S=0;
quiver(X(1:under_samp:end)*1e3,Z(1:under_samp:end)*1e3,vx_LS_Bone(1:under_samp:end,1:under_samp:end,time_index).*scaling_factor*1e3.*PD_threshold_mask(1:under_samp:end,1:under_samp:end,time_index),vz_LS_Bone(1:under_samp:end,1:under_samp:end,time_index).*scaling_factor*1e3.*PD_threshold_mask(1:under_samp:end,1:under_samp:end,time_index),S,'LineWidth',1,'Color','b','MaxHeadSize',0.99)
hold on
plot(X*1e3,fit_curve_Periosteum*1e3,'c-.','linewidth',1)
plot(X*1e3,fit_curve_Endosteum*1e3,'c-.','linewidth',1)
hold off
axis image 
title('LS Vector Doppler')
xlim([X(1) X(end)]*1e3)
ylim([Z(1) Z(end)]*1e3)
axis ij
colormap hot
    xlabel('width [mm]')
    ylabel('depth [mm]')
set(gca,'Fontsize',12)
CurrentPosition = get(gca,'Position'); %[left bottom width height]
NewPosition = CurrentPosition; %[left bottom width height]
NewPosition(2) = NewPosition(2)*1.1;
set(gca,'Position',NewPosition)

figure(111)
subplot(5,1,[3:4])
s2=imagesc(X*1e3,Z*1e3,PowerDoppler_dB(:,:,time_index));
caxis([PD_dB_Threshold 0])
hold on
S=0;
quiver(X(1:under_samp:end)*1e3,Z(1:under_samp:end)*1e3,vx_TO_Bone(1:under_samp:end,1:under_samp:end,time_index).*scaling_factor*1e3.*PD_threshold_mask(1:under_samp:end,1:under_samp:end,time_index),vz_TO_Bone(1:under_samp:end,1:under_samp:end,time_index).*scaling_factor*1e3.*PD_threshold_mask(1:under_samp:end,1:under_samp:end,time_index),S,'LineWidth',1,'Color','b','MaxHeadSize',0.99)
hold on
plot(X*1e3,fit_curve_Periosteum*1e3,'c-.','linewidth',1)
plot(X*1e3,fit_curve_Endosteum*1e3,'c-.','linewidth',1)
hold off
axis image 
title('TO Vector Doppler')
xlim([X(1) X(end)]*1e3)
ylim([Z(1) Z(end)]*1e3)
axis ij
colormap hot
    xlabel('width [mm]')
    ylabel('depth [mm]')
set(gca,'Fontsize',12)
CurrentPosition = get(gca,'Position'); %[left bottom width height]
NewPosition = CurrentPosition; %[left bottom width height]
NewPosition(2) = NewPosition(2)*1.15;
set(gca,'Position',NewPosition)

subplot 515
plot(time,PD,'r')
hold on
plot(time(time_index-index_range_display(1)+1),PD(time_index-index_range_display(1)+1),'k.','markersize',24)
hold off
ylabel({'blood volume','in motion [a.u.]'})
xlabel('time [s]')
xlim([time(1) time(end)])
set(gca,'Fontsize',12)


set(gcf,'Color','w')

pause(0.01)
F(t_ind)=getframe(H1);

t_ind = t_ind +1;

end


%%

% v = VideoWriter('TEST.avi');
v = VideoWriter('TEST.mp4','MPEG-4');
v.Quality = 100;
v.FrameRate = 15;
open(v)
writeVideo(v,F)
close(v)




