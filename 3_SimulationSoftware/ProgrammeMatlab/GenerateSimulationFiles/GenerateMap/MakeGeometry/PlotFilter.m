std = 3;
nhood = 5;
map_gauss_3 = imgaussfilt(uint8(Map), std);
map_bw = bwareaopen(map_gauss_3,1000);
X = 0:0.01:30-0.01; X = X -mean(X);
Z = flipud(0:0.01:15-0.01);
% MAP T Original 
figure, 
% subplot(3, 1, 1) 
imagesc(X(1000:2000), Z(800:1100), Map(800:1100, 1000:2000));
title('\textbf{Simulation map}', 'Generated map', 'Interpreter', 'latex',  'FontSize', 22);

ax = gca;
set(ax, 'Units', 'normalized'); 
ax.FontSize = 16;
axis image
% Gaussian, std = 3
% subplot(3, 1, 2) 
figure,
imagesc(X(1000:2000), Z(800:1100), map_gauss_3(800:1100, 1000:2000));
title('\textbf{Simulation map}', 'Generated map filtred with a Gaussian filter', 'Interpreter', 'latex',  'FontSize', 22);

ax = gca;
set(ax, 'Units', 'normalized'); 
ax.FontSize = 16;
axis image

% Gaussian, std = 2 + voisins 1000
% subplot(3, 1, 3) 
figure,
imagesc(X(1000:2000), Z(800:1100), map_bw(800:1100, 1000:2000));
title('\textbf{Simulation map}', 'Removal of unconnected component + Gaussian filter', 'Interpreter', 'latex',  'FontSize', 22);

ax = gca;
set(ax, 'Units', 'normalized'); 
ax.FontSize = 16;
axis image
%%
[Rq, Corr, rugosity] = ComputeInterfaceParameters(map_bw, grid, probe, medium);

format = '\nThe interface have a rms of %.3f mm and a correlation length of %.3f mm \nThe rugosity in one wavelength is %.1f%%';
fprintf(format,Rq, Corr, rugosity);
