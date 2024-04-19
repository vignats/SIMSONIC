Map;
std = 10;
nhood = 1000;
map_gauss_3 = imgaussfilt(uint8(Map), std);
map_bw = bwareaopen(map_gauss_3,1000);

% MAP T Original 
figure, 
subplot(3, 1, 1) 
imagesc(Map(800:1100, 1000:2000));

title('Initial');


% Gaussian, std = 3
subplot(3, 1, 2) 
imagesc(map_gauss_3(800:1100, 1000:2000));
title(['Gaussian, std =', num2str(std)]);


% Gaussian, std = 2 + voisins 1000
subplot(3, 1, 3) 
imagesc(map_bw(800:1100, 1000:2000));
title(['Gaussian, std =', num2str(std), ' and nhood = ', num2str(nhood)]);

%%
[Rq, Corr, rugosity] = ComputeInterfaceParameters(map_bw, grid, probe, medium);

format = '\nThe interface have a rms of %.3f mm and a correlation length of %.3f mm \nThe rugosity in one wavelength is %.1f%%';
fprintf(format,Rq, Corr, rugosity);
