function [corrRecovered] = ComputeCorr(heights, gridStep)
% This function allows to compute the correlation length of a profile.
% 
% INPUT : heights - vector column that contains the heigth regarding a
% surface planar at origin.
%         step - spatial spacing between consecutive points. (usually
%         grid.step)
% OUTPUT : correlation length of the profile 
% 
% see also : autocorr()

    % Computation of the correlation and lags.
    [acf,lags] = autocorr(heights, 'NumLags', length(heights)-1);

    % The correlation length is the lag that corresponds to a correlation of 1/e
    % Find the closest value from 1/e.
    diffValues = abs(acf - 1/exp(1));
    lagsCorrLength = lags(diffValues == min(diffValues));

    % To obtain the value in mm, we multiply by the heights profil step
    corrRecovered = lagsCorrLength*gridStep;
end