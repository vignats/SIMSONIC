function [f,x] = rsgene1D(N,rL,h,cl)
%
% [f,x] = rsgene1D(N,rL,h,cl) 
%
% generates a 1-dimensional random rough surface f(x) with N surface points.
% The surface has a Gaussian height distribution and an
% exponential autocovariance function, where rL is the length of the surface, 
% h is the RMS height and cl is the correlation length.
%
% Input:    N   - number of surface points
%           rL  - length of surface
%           h   - rms height
%           cl  - correlation length
%
% Output:   f   - surface heights
%           x   - surface points
%
% Last updated: 2010-07-26 (David Bergström).  
%

format long;

x = linspace(-rL/2,rL/2,N);

Z = h.*randn(1,N); % uncorrelated Gaussian random rough surface distribution
                   % with rms height h
                        
% Gaussian filter
F = exp(-abs(x)/(cl/2));

% correlated surface generation including convolution (faltning) and inverse
% Fourier transform and normalizing prefactors
f = sqrt(2)*sqrt(rL/N/cl)*ifft(fft(Z).*fft(F));