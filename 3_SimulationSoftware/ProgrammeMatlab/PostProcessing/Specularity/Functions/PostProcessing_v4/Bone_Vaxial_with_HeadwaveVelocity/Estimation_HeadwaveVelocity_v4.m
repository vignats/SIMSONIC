clear all
close all
clc

%%-----------------------------------------------------------------------%%
%%---- WARNING: before running this code, you must compile the code to
%%---- calculate travel time of specular reflection at periosteum,
%%---- the compiled MEX file is called specular_reflection_PP
%%-----------------------------------------------------------------------%%

displayed_dynamic_range_dB = 50;
filtering_RF_data_YES = 1;
apodization_start_end_YES = 1;
TGC_YES = 1;
correction_delay_matching_layers_YES = 1;

%----------
%----------

acquisition_name = 'Longi_Sujet4_Tibia_N1_HD_1.h5';

% wavespeed model (Here we only reconstruct a soft-tissue image)
C_TISSUE            = 1600; % [m/s]
C_AXIAL             = C_TISSUE; % [m/s]
C_RADIAL            = C_TISSUE; % [m/s]
ANISO_SHAPE_COEF    = 1.0; % unitless
C_MARROW            = C_TISSUE; % [m/s]
%----------
%----------


%%
%%---------------------------------------------------------------------------------%%
%%---- RECON WITH SINGLE-ELEMENT TRANSMISSIONS FOR SEGMENTATION OF PERIOSTEUM  ----%%
%%---------------------------------------------------------------------------------%%


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

SIG = SIG(1:end/2,:,:);
LENGTH = size(SIG,1);

NCORE                   = 8; % number of CPU cores for OpenMP parallel computing

%--------------------------------------------------%
%           P4-1 probe
%--------------------------------------------------%
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
zmin    = 2.5e-3; % [m]
zmax    = 10e-3; % [m]
pixel_size = C_TISSUE/(FREQ_Transducteur*12); % [m]

% create image axes/pixels coordinates
X       = xmin:pixel_size:xmax; % image width [m]
Z       = zmin:pixel_size:zmax; % image depth [m]

% Sources to be used during reconstruction
Tx_to_be_used       = 1:NSOURCES;

% parameter only for ImageTissue, ImageBone, ImageMarrow:
HalfOpeningAngInSkinDeg = 40; % half opening angle in receive only for image reconstructed with full aperture in receive [degree]
HalfOpeningAngInLensRad = asin(C_LENS/C_TISSUE*sind(HalfOpeningAngInSkinDeg)); % [rad]
SubApertureApodis       = 1; % 0, 1, 2 for no windowing, hamming, hann.



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
   center_freq_bandwidth = 2.4e6;
b = fir1(60,[0.6 1.4]*center_freq_bandwidth/(Fs/2),'bandpass'); % FullBand
% b = fir1(110,[0.82 1.18]*center_freq_bandwidth/(Fs/2),'bandpass'); % NarrowBand
    for ii=1:size(SIG,2)
        for jj=1:size(SIG,3)
        SIG(:,ii,jj) = filtfilt(b,1,SIG(:,ii,jj));
        end
    end
    RF_filter = 'fir1;filtfilt;order=110;Wn=[0.6 1.4]';
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

AccuracyBrent = 5e-4;
IterMaxBrent = 30;
nTestPriorTx = 30;
nTestPriorRx = 30;
grid_size_TravelTime = 0;%pixel_size*4; % (m) grid size for interpolation of travel times and reconstructed image
grid_size_Image = 0;%pixel_size*2; % (m)

%-----------%
%-- FLOAT --%
%-----------%
XR = single(XR);
ZR = single(ZR);
XS = single(XS);
ZS = single(ZS);
X = single(X);
Z = single(Z);
I_SIG = single(I_SIG);
Q_SIG = single(Q_SIG);
Tx_to_be_used = single(Tx_to_be_used);
add_to_delay_firing = single(add_to_delay_firing);
%-----------%
%-----------%
        
Setup = [LENS_THICKNESS D_PROBE_PERIOS MIN_CORTICAL MAX_CORTICAL...
    HalfOpeningAngInLensRad Fs C_LENS C_TISSUE C_MARROW C_AXIAL C_RADIAL...
    ANISO_SHAPE_COEF NCORE grid_size_TravelTime grid_size_Image...
    AccuracyBrent IterMaxBrent nTestPriorTx nTestPriorRx SubApertureApodis...
    XMinSegmentationEndo XMaxSegmentationEndo];
Setup = single(Setup);

[FullImage, PeriParab, EndoParab] =...
        ReconTBM_HYPERFAST_FLOAT_Stitch_RxApod_SegEndo_v4(Setup, XR, ZR, XS, ZS, X, Z,...
        I_SIG, Q_SIG, Tx_to_be_used, add_to_delay_firing);
        

Im = 20*log10(FullImage);
Im = Im - max(Im(:));

fit_Periosteum  = PeriParab(1)*X.^2 + PeriParab(2)*X + PeriParab(3);

    
%%

    
rect = [50 50 1000 500]; %[left bottom width height];
H1=figure(88);
set(88,'Position',rect)


%------------------------------------------------------------%
%------------------- HEADWAVE VELOCITY ----------------------%
%------------------------------------------------------------%

C_AXIAL_good_guess = 3900;


% distance from point (x0,y0) to line ax+by+c=0
% dist = abs(a*x0+b*y0+c)/sqrt(a^2+b^2);
a = PeriParab(2);
b = -1;
c = PeriParab(3); 
x0 = XR;
y0 = 0;

dist = abs(a*x0+b*y0+c)/sqrt(a^2+b^2);
n = [a b]/sqrt(a^2+b^2);



dip_angle = atan(PeriParab(2));
margin = 1/FREQ_Transducteur*2;

nb_pts_fft = 2^nextpow2(Nt*2);
frek = [[0:nb_pts_fft/2-1] [-nb_pts_fft/2:-1]]*Fs/nb_pts_fft;
frek = repmat(frek.',1,NELEMENTS);

    
mask_muting = zeros(Nt,size(SIG,2));
count = 1;
count_display = 1;
c_min = C_AXIAL_good_guess*0.6; % m/s
c_max = C_AXIAL_good_guess*1.4; % m/s
wavespeed = c_min:20:c_max;
S1 = zeros(length(wavespeed),2);

for el_ind = [1 NELEMENTS]

    im = squeeze(SIG(:,:,el_ind));
        
    X0 = XS(el_ind);
    ReceiverXPos = XS - X0; % Make transmit element zero

    t_tissue = 2*(dist(el_ind)-LENS_THICKNESS/abs(n(2)))/C_TISSUE;

    ic_tissue_bone = asin(C_TISSUE/C_AXIAL_good_guess);
    temp1 = asin(C_LENS/C_TISSUE*sin(ic_tissue_bone-dip_angle));
    temp2 = asin(C_LENS/C_TISSUE*sin(ic_tissue_bone+dip_angle));    
    arrival_time_headwave = t_tissue*cos(ic_tissue_bone) + LENS_THICKNESS/C_LENS*(1/cos(temp1)+1/cos(temp2)) + (abs(ReceiverXPos)-LENS_THICKNESS*(tan(temp1)+tan(temp2)))/C_TISSUE.*(sin(ic_tissue_bone+dip_angle)*(XS>X0) + sin(ic_tissue_bone-dip_angle)*(XS<X0));
    arrival_time_headwave = arrival_time_headwave.*(abs(ReceiverXPos) > 4.5e-3);
    
    C = double([C_LENS C_TISSUE]);
    P = double([PeriParab(1:3) LENS_THICKNESS]);
    XR = double(XR);
    XS = double(XS);
    [travel_timePP, ~, ~] = specular_reflection_PP(C, P, XS(el_ind), XR, 0, 1);    
    
    
    for jj=1:NELEMENTS
        mask_muting(:,jj) =  (Time  < travel_timePP(jj) - 1*margin) & (Time > arrival_time_headwave(jj)-2.*margin) & (Time < arrival_time_headwave(jj)+2.*margin) & ((travel_timePP(jj) - arrival_time_headwave(jj)) > 1.0*margin) & (Time > 4*LENS_THICKNESS/C_LENS);%& (abs(XS(el_ind)-XS(jj))>10e-3);%
    end
        
        
    figure(88)
    subplot(2,4,count_display)
    imagesc(XR*1e3,Time*1e6,im)
    hold on
    plot(XS(abs(ReceiverXPos) > 4e-3)*1e3,arrival_time_headwave(abs(ReceiverXPos) > 4e-3)*1e6,'r-.','linewidth',2)
    plot(XS*1e3,travel_timePP*1e6,'y-','linewidth',2)
    hold off
    caxis([-1 1]*max(abs(im(:)))/10)
    xlabel('azimutal distance (mm)')
    ylabel('arrival time (\mus)')
    title(['Tx with element ' int2str(el_ind)])
    set(gca,'fontsize',12)
    colormap gray
            
    im = im.*mask_muting;
    im = im./repmat(max(max(im),ones(1,NELEMENTS)),Nt,1);
    
    figure(88)
    subplot(2,4,count_display+4)
    imagesc(XR*1e3,Time*1e6,im)
    hold on
    plot(XS(abs(ReceiverXPos) > 4e-3)*1e3,arrival_time_headwave(abs(ReceiverXPos) > 4e-3)*1e6,'r-.','linewidth',2)
    plot(XS*1e3,travel_timePP*1e6,'y-','linewidth',2)
    hold off
    caxis([-1 1]*max(abs(im(:)))/5)
    xlabel('azimutal distance (mm)')
    ylabel('arrival time (\mus)')
    title('muted echo signals')
    set(gca,'fontsize',12)



        
    %---------- Semblance analysis --------------%

    for ii=1:length(wavespeed)

        Tx_delay_per_pitch = PITCH/wavespeed(ii);
        delay_time = [0:NELEMENTS-1]*Tx_delay_per_pitch;
        sig_delayed = zeros(nb_pts_fft,NELEMENTS);
        apod_mat = zeros(nb_pts_fft-Nt,NELEMENTS);
        RF = cat(1,apod_mat,im);
        FFT_matrix = fft(RF,nb_pts_fft);
        if count==1
            sig_delayed = pulse_delayingXFFT(FFT_matrix,-delay_time,frek);
        elseif count==2
            sig_delayed = pulse_delayingXFFT(FFT_matrix,-fliplr(delay_time),frek);
        end
       S1(ii,count) = semblance(sig_delayed);

    end


    [~,index_max] = max(S1(:,count));
    if (index_max==1)||(index_max==length(wavespeed))
        headwave_velocity(count) = wavespeed(index_max);
    else
    peak = S1(index_max-1:index_max+1,count);
    y1 = peak(1); % y = ax2^+ bx + c
    y2 = peak(2);
    y3 = peak(3);
    a=(y1-2*y2+y3)/2.0;
    b=-(3*y1-4*y2+y3)/2.0 - 2.0  *a*(index_max-1);
    c=y1-a*(index_max-1)*(index_max-1)-b*(index_max-1);
    index_fraction = (- b/2/a );
    headwave_velocity(count) = wavespeed(floor(index_fraction)) + (index_fraction-floor(index_fraction))*(wavespeed(2)-wavespeed(1));
    %subsample_ampl = a*(- b/2/a )^2 + b*(- b/2/a ) + c;
    end
   
    count = count+1;
    count_display = count_display+1;
end


    
%c_bone_semblance_small_dip_angle = 2*headwave_velocity(1)*headwave_velocity(2)/(headwave_velocity(1)+headwave_velocity(2)) % valid if dip angle < 6? or 0.1 rad
c_bone_semblance = 2*headwave_velocity(1)*headwave_velocity(2)/(headwave_velocity(1)+headwave_velocity(2))*cos(dip_angle);

figure(88)
subplot(2,4,7:8)
plot(wavespeed,S1)
xlim([c_min c_max])
ylabel('wavefront coherence')
legend('Tx with element 1',['Tx with element ' int2str(NELEMENTS)])
xlabel('tested wavespeed values [m/s]')
title(['headwave velocity = ' num2str(c_bone_semblance,4) ' m/s']);
set(gca,'fontsize',12)

figure(88)
subplot(2,4,3:4)
imagesc(X*1e3,Z*1e3,Im)
ylabel('Depth [mm]')
xlabel('Width [mm]')
hold on
plot(X*1e3,fit_Periosteum*1e3,'r')
hold off
axis image
caxis([-30 0]);
title('Soft tissue image')
set(gca,'fontsize',12)

