function [VZ_Full,VZ_TO,VX_TO,RFallHFR1,RFallHFR2]=TO_Dec2022(ImIQ,param)


sigmax=param.sigmax; % [mm]
X=param.X; % [m]
Z=param.Z; % [m]
dx=(X(2)-X(1))*1000; % [mm]
dz=(Z(2)-Z(1))*1000; % [mm]
lag=param.lag;

Nz=size(ImIQ,1);
Nx=size(ImIQ,2);
Nt=size(ImIQ,3);



f0x = 1/param.lambdax; % [1/mm]
fx = linspace(-1/(2*dx),+1/(2*dx),Nx);
% fz = linspace(-1/(2*dz),+1/(2*dz),Nz);

% create gaussian mask
rightAperture=exp(-2*(sigmax*pi*(fx-f0x)).^2);
leftAperture =fliplr(rightAperture);

leftAperture_mask=repmat(leftAperture,[Nz,1,Nt]);  %kspace aperture for whole dataset
rightAperture_mask=repmat(rightAperture,[Nz,1,Nt]);

% FFT of data
IQimFFT=fft((ImIQ),[],2); % fft in lateral dimension


%---- TO
% apply left-lobe apodization to data (freq domain)
RFallHFR1=ifft(IQimFFT.*ifftshift(leftAperture_mask),[],2);
% apply right-lobe apodization to data (freq domain and back)
RFallHFR2=ifft(IQimFFT.*ifftshift(rightAperture_mask),[],2);
% % add apertures together
% RFOT=RFallHFR1+RFallHFR2;




% Sampling
Nkx = 2^nextpow2(Nx)*4;
Nkz = 2^nextpow2(Nz)*4;

% Frequency and wave number vectors
fx = [[-Nkx/2:-1] [0:Nkx/2-1]]/(dx)/Nkx;
fz = [[-Nkz/2:-1] [0:Nkz/2-1]]/(dz)/Nkz;

%---- estimate actual transverse spatial frequency in Fourier domain
tmp = abs(fftshift(fft2(squeeze(RFallHFR1(:,:,1)),Nkz,Nkx)));
% tmp = abs(IQimFFT.*ifftshift(leftAperture_mask));
fx_grid = repmat(fx,[Nkz,1]);
tmp2 = fx_grid.*tmp.^2;
actual_f0x_leftAperture = sum(tmp2(:))/sum(tmp(:).^2);
tmp = abs(fftshift(fft2(squeeze(RFallHFR2(:,:,1)),Nkz,Nkx)));
% tmp = abs(IQimFFT.*ifftshift(rightAperture_mask));
% fx_grid = repmat(fx,[Nz,1,Nt]);
tmp2 = fx_grid.*tmp.^2;
actual_f0x_rightAperture = sum(tmp2(:))/sum(tmp(:).^2);

disp(['estimated lateral spatial period in Fourier (left aperture)=' num2str(1/actual_f0x_leftAperture,3) ' mm'])
disp(['estimated lateral spatial period in Fourier (right aperture)=' num2str(1/actual_f0x_rightAperture,3) ' mm'])



if(param.display)
figure(12);
xx=[0:size(ImIQ,2)-1]*dx;
zz=[0:size(ImIQ,1)-1]*dz;
xx = xx - mean(xx);
subplot(221);
%imagesc(xx,zz,real(ImIQ(:,:,1)));
imagesc(xx,zz,real(RFallHFR1(:,:,1)));
axis image
%title('Initial Data')
title('left aperture')
xlabel('[mm]')
ylabel('[mm]')
subplot(222);
%imagesc(xx,zz,real(RFOT(:,:,1)));
imagesc(xx,zz,real(RFallHFR2(:,:,1)));
axis image
%title('TO Data')
title('right aperture')
xlabel('[mm]')
ylabel('[mm]')


subplot(223);
im = 20*log10(abs(fftshift(fft2(squeeze(ImIQ(:,:,1)),Nkz,Nkx))));
im = im - max(im(:));
imagesc(fx,fz,im);
hold on
plot(f0x*[1 1],[fz(1) fz(end)],'r')
plot(-f0x*[1 1],[fz(1) fz(end)],'r')
hold off
caxis([-50 0])
axis image
colorbar
xlabel('kx [1/mm]')
ylabel('kz [1/mm]')
title('Initial Spectrum')
subplot(224);
im = 20*log10(abs(fftshift(fft2(squeeze(RFallHFR1(:,:,1)+(RFallHFR2(:,:,1))),Nkz,Nkx))));
im = im - max(im(:));
imagesc(fx,fz,im);
hold on
plot(f0x*[1 1],[fz(1) fz(end)],'r')
plot(-f0x*[1 1],[fz(1) fz(end)],'r')
hold off
caxis([-50 0])
axis image
colorbar
xlabel('kx [1/mm]')
ylabel('kz [1/mm]')
title('TO Spectrum')

%pause
end


%%Autorrelation
% autocorrelate data one frame with the next
RZ=ImIQ;
RZ=RZ(:,:,1:end-lag).*conj(RZ(:,:,lag+1:end))/lag;

% autocorrelate left apodized data one frame with the next
R1 = (RFallHFR2);
R1 = R1(:,:,1:end-lag).*conj(R1(:,:,lag+1:end))/lag;

% autocorrelate right apodized data one frame with the next
R2 = (RFallHFR1);
R2 = R2(:,:,1:end-lag).*conj(R2(:,:,lag+1:end))/lag;

if (param.Nt~=1)&&(param.Nt~=0)
    %-----------------------------------------%
    %------- temporal smoothing -----%
    %-----------------------------------------%
    NT = param.Nt;
    ht=ones(1,NT)/NT;
    RZ = filter(ht,1,RZ,[],3);
    R1 = filter(ht,1,R1,[],3);
    R2 = filter(ht,1,R2,[],3);
end

if (param.Nt==0)
    %----------------------------------------------%
    %------- temporal averaging on all frames -----%
    %----------------------------------------------%
    RZ = mean(RZ,3);
    R1 = mean(R1,3);
    R2 = mean(R2,3);
end

if (param.Nx~=1 && param.Nz~=1)
    %-----------------------------------------%
    %------- spatial averaging/smoothing -----%
    %-----------------------------------------%
    Nx = param.Nx;
    Nz = param.Nz;
    sAvg = ones(Nz,Nx)./Nx/Nz;
    for i=1:size(RZ,3)    
        RZ(:,:,i)= filter2(sAvg,RZ(:,:,i)); 
        R1(:,:,i)= filter2(sAvg,R1(:,:,i)); 
        R2(:,:,i)= filter2(sAvg,R2(:,:,i)); 
    end
end


[~,median_vertical_spatial_period_image,~,median_lateral_spatial_period_image,spatial_period]=estimate_spatial_period_image(ImIQ,X,Z,param.assumed_spatial_period_image,param.mask_ROI,param.display);
median_vertical_spatial_period_image = median_vertical_spatial_period_image*1000; % in mm
disp(['estimated vertical spatial period in image (without TO)=' num2str(median_vertical_spatial_period_image,3) ' mm'])
% compute VZ with 0 degree PSF direction (traditional axial Doppler)
VZ_Full=angle(RZ)/(2*pi)*median_vertical_spatial_period_image*param.fps; % same as v1*(1/(4*pi*param.f0))*param.fps*param.c0;

% compute VZ and VX with TO
v1=angle(R1);
v2=angle(R2);

vxTO=(v1-v2)/2;
vzTO=(v1+v2)/2;


[~,median_vertical_spatial_period_image,~,median_lateral_spatial_period_image,spatial_period]=estimate_spatial_period_image(RFallHFR1,X,Z,param.assumed_spatial_period_image,param.mask_ROI,param.display);
median_vertical_spatial_period_image = median_vertical_spatial_period_image*1000; % in mm
median_lateral_spatial_period_image = median_lateral_spatial_period_image*1000; % in mm
disp(['estimated vertical spatial period in image (with TO)=' num2str(median_vertical_spatial_period_image,3) ' mm'])
disp(['estimated lateral spatial period in image (with TO)=' num2str(median_lateral_spatial_period_image,3) ' mm'])

VZ_TO=vzTO/(2*pi)*median_vertical_spatial_period_image*param.fps; % mm/s
VX_TO=vxTO/(2*pi)*median_lateral_spatial_period_image*param.fps; % mm/s

end