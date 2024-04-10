% This script aim to study the specular probability in a zone of interest
% in order to quantify the response of the specular model for various RMS
% and correlation length. 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% Compute results 
simuDirAll = '/calculSSD/salome/Simulation-28mars';
simuNameW = '/Bone227G-Image1590-F0.06/';
simuNameR = '/Bone227G-Image1590-F1.25/';

[~, ProbaR, OrientationtR, reconstructionR] = PostProcessing([simuDirAll, simuNameR]);
parametersR = load(fullfile(simuDirAll, simuNameR, 'parameters.mat'));

[~, ProbaW, OrientationtW, reconstructionW] = PostProcessing([simuDirAll, simuNameW]);
parametersW = load(fullfile(simuDirAll, simuNameW, 'parameters.mat'));


%% Extract the zone of interest on the specular map probability
plotLimite = true;
nbWavelength = 1;
[ROI_R, LimInfR, LimSupR] = ExtractROI(ProbaR, parametersR, reconstructionR, nbWavelength, plotLimite);
[ROI_W, LimInfW, LimSupW] = ExtractROI(ProbaW, parametersW, reconstructionW, nbWavelength, plotLimite);

%% Plot ROI 
figure
subplot(2,1,1)
pcolor(ROI_R) 
shading flat
xlabel('Lateral position (mm)', Interpreter='latex')
ylabel('Depth (mm)', Interpreter='latex')
title('Probability map', 'Roughness profil') 
axis ij image;
colorbar
subplot(2,1,2)
pcolor(ROI_W) 
shading flat
xlabel('Lateral position (mm)', Interpreter='latex')
ylabel('Depth (mm)', Interpreter='latex')
title('Probability map', 'Waviness profil') 
axis ij image;
colorbar

%% COMPUTE METRICS 
linearProbaR = mean(ROI_R, 1);
totalProbaR = mean(linearProbaR);

linearProbaW = mean(ROI_W, 1);
totalProbaW = mean(linearProbaW);

figure;
plot(reconstructionR.Xmm, linearProbaR, 'b');
hold on 
plot(reconstructionR.Xmm, linearProbaW, 'r');
xlabel('Lateral position (mm)', Interpreter='latex')
ylabel('Specular probability', Interpreter='latex')
title('Specular probability along the lateral position', sprintf('Mean probability along roughness profile is %.02f and waviness profile is %.02f', totalProbaR, totalProbaW));
legend('Roughness', 'Waviness')
ylim([0, 1]);
%%
function[ROI, LimInf, LimSup] = ExtractROI(SpecularProbaMap, parameters, reconstruction, nbWavelength, plotLimit)
% This function allows to extract the zone of interest that coresponds to
% the bone/tissu interface in the specularity probability map.

    % Interface depth
    to_px = @(mm) round(mm*1e-3/reconstruction.pixel_size)+1;       % The values in reconstruction are expressed in meter
    specularDepth = 1.7;                                              % The specularity Map starts at 2mm from the probe.
    specularInterface = parameters.interface.depth - parameters.probe.depth - specularDepth;
    
    % Length of the ROI
    lengthROI = nbWavelength*(parameters.medium.cp(2) / (parameters.probe.fc)*1e3);       % Number of wavelength in the bone that we consider (mm)
    LimInf = to_px(specularInterface - lengthROI/2);
    LimSup = to_px(specularInterface + lengthROI/2);
    ROI = SpecularProbaMap(LimInf : LimSup, :);
    
    if plotLimit
        figure
        pcolor(reconstruction.Xmm, reconstruction.Zmm, SpecularProbaMap) 
        yline(reconstruction.Zmm(LimInf), 'Color', 'red', 'LineWidth', 2);
        yline(reconstruction.Zmm(LimSup) , 'Color', 'red', 'LineWidth', 2);
        shading flat
        xlabel('Lateral position (mm)', Interpreter='latex')
        ylabel('Depth (mm)', Interpreter='latex')
        title('Probability map') 
        axis ij image;
        colorbar
    end
end
