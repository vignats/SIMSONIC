function[SpecularModel, SpecularProbaMap, OrientationtMap, reconstruction] = getPostProcess(simuDir)
% This function returns the specularity model, the probability of
% specularity map and the tilt angle map of the simulation contained in
% simuDir
% 
% Input : simuDir - directory containing the simulations (tx_files)
% Output : SpecularModel - model of perfectly specular signal receive from each pixel
%          SpecularProbaMap - map of specular index
%          OrientationtMap - map of the orientations angles
%          reconstruction - structur of the parameters used for the
%          reconstruction.
%
% See also : LoadRfData, GenerateParamRecon,ComputeSpecularTransform, ComputeSpecularityModel

    % Get parameters of the simulation 
    parameters = load(fullfile(simuDir, 'parameters.mat'));
    recorded = LoadRfData(parameters.probe, simuDir);
    
    % Compute parameters required to reconstruct the image using DAS and/or specular transform
    [acquisition, reconstruction] = GenerateParamRecon(recorded, parameters, simuDir);
    
    % Compute the specular transform
    [SpecularTransform, TiltAngles] = ComputeSpecularTransform(reconstruction, acquisition, parameters);
    
    % Compute teh specularity map and angle
    plotMap = false;
    [SpecularModel, SpecularProbaMap, OrientationtMap] = ComputeSpecularityModel(SpecularTransform, ...
        acquisition, reconstruction, TiltAngles, plotMap, parameters, simuDir);
end