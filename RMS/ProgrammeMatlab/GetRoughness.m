function[profile, roughness] = GetRoughness(filename, fc, segmented)
% This function allow to filter the low frequency data of the profile,
% which corresponds to the waviness, in order to obtain the roughness that
% corresponds to higher frequency. 
% INPUT :   filename - path to the image.
%           fc - Cutt-off frequency between waviness and roughness (mm-1).
%           segmented - indicate weither the image is already segmented or
%           not. 
% OUTPUT : roughness - 1-D vector containing the heigth of the roughness (mm).
%
% See also : ExtractBoundary
    
    % BOUNDARIES COMPUTATION
    boundary_endost = ExtractBoundary(filename,segmented);
    boundary_endost(2, :) = max(boundary_endost(2, :)) - boundary_endost(2, :);
    
    % ROUGHNESS COMPUTATION
    % Parameters
    delta_x = 9e-3;             % Pixel size corresponding to the space step (mm)
    Fs = 1/delta_x;             % Sampling frequency (mm-1)
    profile = boundary_endost(2, :) * delta_x;
    
    roughness = highpass(profile, fc, Fs);
end