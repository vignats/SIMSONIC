addpath(genpath('~/Documents'));

% Pathway to the bone image to process
bones = {'245D', '227G', '267G'};
slices = {{1134, 3195, 3852, 5511}, {2002, 3595, 5614, 6721}, {1700, 3410, 5478, 6716}}; 

b = 3;
s = 1;
dirName = '~/Documents/BoneRugosity/2_ParametersDefinition/BoneImage';
fileName = sprintf('SAMPLE_%s_SLICE_%04d.bmp', bones{b}, slices{s}); 
filePath = fullfile(dirName, fileName);
% Load image
bone_bmp = imread(filePath); 

%% INITIAL COMPUTATION OF RMS HEIGTH
% THRESHOLD TO CONVERT TO BINARY IMAGE
threshold = graythresh(bone_bmp); % Find an automatic threshold
binary_image = imbinarize(bone_bmp, threshold);

% BOUNDARIES COMPUTATION
boundary_endost = ExtractBoundary(binary_image);

figure
imshow(binary_image);
hold on;
plot(boundary_endost(1, :), boundary_endost(2, :), '*');  % Trace le contour
hold on;
xlabel('Width');
ylabel('Depth');
title('Endost boundary'); hold off;

boundary_endost(2, :) = max(boundary_endost(2, :)) - boundary_endost(2, :);
%% COMPUTE FOURRIER TRANSFORM OF THE PROFILE
% Parameters
delta_x = 9e-3;             % Pixel size corresponding to the space step (mm)
Fs = 1/delta_x;             % Sampling frequency (mm-1)
L = size(boundary_endost, 2);
% L = 2048;
% Initial profile
profile = boundary_endost(2, :) * delta_x; % - mean(boundary_endost(2, :));

x = (0:L-1)*1/Fs;

% figure
% plot(x, profile)
% title('Original profile with rugosity');
% xlabel('Width (mm)');
% ylabel('Depth (mm)');

% Fourrier transform
spectrum = fft(profile);

% f = Fs/L*(0:L-1);
f = Fs/L*(-L/2:L/2 - 1);
figure
plot(f, fftshift(abs(spectrum))) 
title("Spectrum")
xlabel("f (mm^{-1})")
ylabel("|fft(X)|")

scaled_spectrum = abs(spectrum/L);
scaled_spectrum = scaled_spectrum(1:floor(L/2) + 1);
scaled_spectrum(2:end-1) = 2*scaled_spectrum(2:end-1);

f_side = Fs/L*(0:(floor(L/2)));
figure
plot(f_side,scaled_spectrum,"LineWidth",3) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|FFT(f)|")

%% POWER SPECTRAL DENSITY
figure
plot(f_side,pow2db(scaled_spectrum))
grid on
title("Periodogram Using FFT")
xlabel("Frequency (mm-1)")
ylabel("Power/Frequency (dB*mm)")

%% FILTER THE HIGH FREQUENCY DATA
filtered_spectrum = spectrum(abs(spectrum) > 1);
% filtered_spectrum = spectrum(1:100); 
filtered_profile = ifft(filtered_spectrum);

% Interpolate filtered profile to match length
filtered_profile = interp1(linspace(min(x), max(x), numel(filtered_profile)), filtered_profile, x, 'linear');

% Plot normalized profiles
figure
plot(x, real(filtered_profile)/max(real(filtered_profile)), 'b', 'DisplayName', 'Filtered Profile');
hold on
% figure
plot(x, profile/max(profile), 'r', 'DisplayName', 'Original Profile');
title('Comparison of Original and Filtered Profiles');
xlabel('Width (mm)');
ylabel('Depth (Normalized)');
legend();

%% WITH BUILT-IN FUNCTION 
fpass = 1.25;
waviness = lowpass(profile,fpass,Fs);

figure
plot(x, real(waviness)/max(real(waviness)), 'b', 'DisplayName', 'Waviness');
hold on
plot(x, profile/max(profile), 'r', 'DisplayName', 'Profile');
title('Comparison of Original and Filtered Profiles');
xlabel('Width (mm)');
ylabel('Depth (Normalized)');
legend();

roughness = highpass(profile, fpass, Fs);
figure
plot(x, roughness, 'r');
title('Roughness of the interface');
xlabel('Width (mm)');
ylabel('Depth (mm)');
axis equal

%% FILTERING WITH A MOVING MEAN
profile_mean = movmean(profile, 50);

figure
plot(x, real(profile_mean)/max(real(profile_mean)), 'b', 'DisplayName', 'Moving mean');
hold on
plot(x, profile/max(profile), 'r', 'DisplayName', 'Profile');
title('Comparison of Original and Mean Profiles');
xlabel('Width (mm)');
ylabel('Depth (Normalized)');
legend();
