function [SpecularTransform, TiltAngles] = ComputeSpecularTransform(reconstruction, acquisition, parameters)
    tiltMax = 70; NTilts = 2*tiltMax + 1;
    TiltAngles = linspace(-tiltMax,tiltMax,NTilts);
    fprintf('----------Getting 1D specular transform----------\n');

    if isfield(parameters, 'bone') || (isfield(parameters.interface, 'periost') && parameters.interface.periost > 0)
        fprintf('\tSoft Tissues\n')
        tic
        [SpecularTransform.Tissu] = function_get_specular_transform_interpol(reconstruction, reconstruction.Tissu.timeFlight, ...
            reconstruction.Tissu.Pixel.angleView, reconstruction.Tissu.Degre.angleView, acquisition.Fs, TiltAngles);
        toc
    end
    fprintf('\tBone Tissues\n')
    tic
    % [SpecularTransform.Bone] = function_get_specular_transform_interpol(reconstruction, reconstruction.Bone.timeFlight, ...
    %     reconstruction.Bone.Pixel.angleView, reconstruction.Bone.Degre.angleView, acquisition.Fs, TiltAngles);
    [SpecularTransform.Bone] = get_specular_transform(reconstruction, reconstruction.Bone.timeFlight, ...
        reconstruction.Bone.Degre.angleView, acquisition.Fs, TiltAngles);
    toc
end