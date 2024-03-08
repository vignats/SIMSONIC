clear all, close all , clc


N = 1000 ; %    - number of surface points
rL = 20 ;  % mm - length of surface
h =0.1 ; % mm  - rms height
cl = .5 % mm  - correlation length
% Output:   f   - surface heights
%           x   - surface points
[f,x] = rsgene1D(N,rL,h,cl)

figure, plot(x,f,'.'), axis equal
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