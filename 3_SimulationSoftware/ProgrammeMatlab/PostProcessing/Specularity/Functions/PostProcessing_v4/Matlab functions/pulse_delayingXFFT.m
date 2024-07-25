% correction of the delay in RF line
%
% 3rd May 2010, Version 1
%
% Guillaume Renaud, Erasmus MC, Rotterdam

% a negative delay value moves the signal earlier in time


function [sig_delayed] = pulse_delayingXFFT(FFT_matrix,delay_time,frek)
% first dimension of FFT_matrix is frequency, second is index of channel
% delay_time is a vector

Nb = size(FFT_matrix,1);

delay_time = repmat(delay_time,Nb,1);
% size(FFT_matrix)
% size(frek)
% size(delay_time)
TF_sig_shifted = FFT_matrix.*exp(-1i*2*pi*frek.*delay_time);
sig_delayed = real(ifft(TF_sig_shifted));

return                    
                    
                    
                    
                    










