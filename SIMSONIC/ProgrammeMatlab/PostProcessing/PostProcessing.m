function[SpecularModel, SpecularProbaMap, OrientationtMap, reconstruction] = PostProcessing(simuDir)
% This function returns the specularity model, the probability of
% specularity map and the tilt angle map of the simulation contained in
% simuDir
%
% See also : LoadRfData, GenerateParamRecon,ComputeSpecularTransform, ComputeSpecularityModel

    % Get parameters of the simulation 
    parameters = load(fullfile(simuDir, 'parameters.mat'));
    recorded = LoadRfData(parameters.probe, simuDir);
    
    % Compute parameters required to reconstruct the image using DAS and/or specular transform
    [acquisition, reconstruction] = GenerateParamRecon(recorded);
    
    % Compute the specular transform
    [SpecularTransform, TiltAngles] = ComputeSpecularTransform(reconstruction, acquisition);
    
    % Compute teh specularity map and angle
    plotMap = false;
    [SpecularModel, SpecularProbaMap, OrientationtMap] = ComputeSpecularityModel(SpecularTransform, acquisition, reconstruction, TiltAngles, plotMap, simuDir);
end