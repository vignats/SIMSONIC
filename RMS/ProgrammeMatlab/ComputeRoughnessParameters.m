function [Rq, CorrL] = ComputeRoughnessParameters(roughness)
% ComputeRoughnessParameters compute the root mean square and correlation
% length of a rough surface. 
% INPUT : roughness - high frequency information of a profile. The waviness has been filtered.  
% OUTPUT : Rq - Root mean square of the height roughness profile
%          CorrL - Correlation length of the heigth profile. 

    % Root Mean Square 
    Rq = rms(roughness);

    % Correlation length
    [acf,lags] = autocorr(roughness)
end