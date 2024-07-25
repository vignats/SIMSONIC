function[profile, roughness, xProfile, kc] = FilterProfil(binaryImage, kc)
% This function allow to filter the low frequency data of the profile,
% which corresponds to the waviness, in order to obtain the roughness that
% corresponds to higher frequency. 
% INPUT :   filename - path to the image.
%           fc - Cutt-off frequency between waviness and roughness (mm-1), 
%               if fc = 0, the algorithm find the cutt off frequency.
%           segmented - indicate weither the image is already segmented or
%           not. 
% OUTPUT : roughness - 1-D vector containing the heigth of the roughness (mm).
%
% See also : ExtractBoundary
    
    % BOUNDARIES COMPUTATION
    [boundaryEndost, ~] = ExtractBoundary(binaryImage);
    boundaryEndost(2, :) = max(boundaryEndost(2, :)) - boundaryEndost(2, :);
    
    % ROUGHNESS COMPUTATION
    % Parameters
    delta_x = 9e-3;             % Pixel size corresponding to the space step (mm)
    Fs = 1/delta_x;             % Sampling frequency (mm-1)
    profile = boundaryEndost(2, :) * delta_x;
    xProfile = (boundaryEndost(1,:) - boundaryEndost(1,1))*1/Fs;
    
    if kc == 0
        frequencies = logspace(log10(1e-3), log10(50), 100);
        meanF = zeros(1, length(frequencies));
        
        for i = 1 : length(frequencies)
            filtredProfil = highpass(profile, frequencies(i), Fs);         
            meanF(i) = mean(filtredProfil);
        end
        
        % Find the frequency that corresponds to the roughness of the
        % profile
        kc = frequencies(find(gradient(meanF)./gradient(frequencies) >= 0, 1));
    end
    roughness = highpass(profile, kc, Fs);
end
