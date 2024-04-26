function[probability] = ProbaROI(probaMap, reconstruction, parameters, plotROI)
% This function give the mean specular probability in a zone of interest as
% well as the standard deviation of this probability in order to quantify 
% the response of the specular model for various profile. 
% INPUT : simuDir - directory that contains the simulation of the profil 
% of interest.
%         probaMap - map of the probability of the specularity of the
%         simulated profile
%         reconstruction - parameters of the reconstruction of the image
%         and the computation of the probability of specularity. 
%         parameters - parameters of the simulation.
%         plotROI - boolean that indicate weither to plot the limitation of
%         the ROI. 
%
% OUTPUT : probability - structur that contains the linear probability
% along the lateral profil, the mean probability of specularity over the
% ROI, the standrad deviation of this probability and its correlation length. 
%
% See also : ComputeCorr

    % Extract the region of interest 
    ROI = ExtractROI(probaMap, parameters, reconstruction, plotROI);

    % Compute the metrics of interest
    probability.linearROI = mean(ROI, 1);
    probability.stdROI = std(probability.linearROI);
    probability.meanROI = mean(probability.linearROI);
    probability.corrROI = ComputeCorr(probability.linearROI, reconstruction.pixel_size*1e3);
    probability.stdROIAll = std(ROI);
end

function[ROI] = ExtractROI(SpecularProbaMap, parameters, reconstruction, plotROI)
% This function allows to extract the zone of interest that coresponds to
% the bone/tissu interface in the specularity probability map.

    % Interface depth
    to_px = @(mm) round(mm*1e-3/reconstruction.pixel_size)+1;       % The values in reconstruction are expressed in meter
    specularDepth = 1.7;                                              % The specularity Map starts at 2mm from the probe.
    specularInterface = parameters.interface.endost - parameters.probe.depth - specularDepth;
    
    % Length of the ROI
    lengthROI = (parameters.medium.cp(2) / (parameters.probe.fc)*1e3);       % Number of wavelength in the bone that we consider (mm)
    LimInf = to_px(specularInterface - lengthROI/2);
    LimSup = to_px(specularInterface + lengthROI/2);
    ROI = SpecularProbaMap(LimInf : LimSup, :);
    
    if plotROI
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
