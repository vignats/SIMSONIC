addpath(genpath('~/Documents'));

% Pathway to the bone image to process
bone = {'271G' '267G' '251G' '245D' '227G'};
bone_nb = 5;
dirname = ['~/Documents/BoneRugosity/RMS/', bone{bone_nb}, '/'];
num_slice = 2436;
file = ['SAMPLE_', bone{bone_nb}, '_SLICE_', int2str(num_slice), '.bmp']; 
filename = [dirname , file];

% Load image
bone_bmp = imread(filename); 

% Compute rms and area for different values of cut off frequency
frequencies = logspace(log10(1e-5), log10(55), 100);
parameters = zeros(3, length(frequencies));
for i = 1 : length(frequencies)
    [roughness] = GetRoughness(bone_bmp, frequencies(i));
    parameters(1, i) = rms(roughness);
    parameters(2, i) = trapz(roughness);
    parameters(3, i) = mean(roughness);
end

% Plot evolutions
figure
subplot(1,3,1)
semilogx(frequencies, parameters(1, :));
title('Evolution of the RMS regarding the cutt off frequency');
xlabel('Cutt off frequency (mm-1)');
ylabel('RMS (mm)');

subplot(1,3,2)
semilogx(frequencies, parameters(2, :));
title('Evolution of the area under the boundary regarding the cutt off frequency');
ylabel('Area (mm2)');
xlabel('Cutt off frequency (mm-1)');

subplot(1,3,3)
semilogx(frequencies, parameters(3, :));
title('Evolution of the mean height regarding the cutt off frequency');
ylabel('Mean heigth (mm)');
xlabel('Cutt off frequency (mm-1)');

%% 
fc = [];
for num_slice = randsample((1210 : 2793), 5)
    % Compute rms and area for different values of cut off frequency
    frequencies = logspace(log10(1e-5), log10(55), 100);
    Rq = zeros(1, length(frequencies));
    meanH = zeros(1, length(frequencies));
    for i = 1 : length(frequencies)
        [roughness] = GetRoughness(bone_bmp, frequencies(i));
        Rq(1, i) = rms(roughness);
        meanH(1, i) = mean(roughness);
    end
    
    % Plot evolutions
    figure
    subplot(1,2,1)
    semilogx(frequencies, Rq);
    title('Evolution of the RMS regarding the cutt off frequency', num_slice);
    xlabel('Cutt off frequency (mm-1)');
    ylabel('RMS (mm)');
    
    subplot(1,2,2)
    semilogx(frequencies, meanH);
    title('Evolution of the mean height regarding the cutt off frequency', num_slice');
    ylabel('Mean heigth (mm)');
    xlabel('Cutt off frequency (mm-1)');
end