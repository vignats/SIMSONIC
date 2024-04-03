function [SpecularTransform, TiltAngles] = ComputeSpecularTransform(reconstruction, acquisition)
    tiltMax = 70; NTilts = 2*tiltMax + 1;
    TiltAngles = linspace(-tiltMax,tiltMax,NTilts);
    
    fprintf('----------Getting 1D specular transform----------\n');
    tic
    [SpecularTransform] = function_get_specular_transform_interpol(reconstruction, acquisition, TiltAngles);
    toc
end