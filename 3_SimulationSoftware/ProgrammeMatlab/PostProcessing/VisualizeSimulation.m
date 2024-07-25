% Script to visualize the snapshot of a simulation 
% Parameters computation
simulationPath = '/calculSSD/salome/Simulation-29avr/simulation_rms_0.08_cl_0.5';
tx = 1; 

to_px = @(mm) round(mm/grid.step); 
[Map,Nx,Nz] = SimSonic2DReadMap2D(fullfile(simulationPath, 'Geometry.map2D'));

X = 0:grid.step:grid.width-grid.step; X = X - mean(X);
Z = flipud(0:grid.step:grid.depth-grid.step);
snapNb = param.length/0.5;      % Number of snap per transmitter, one stap recorded every 0.5 us. . 

profile = zeros(1, size(Map, 2));
for i = 1:size(Map, 2)
    profile(1, i) = find(Map(:, i) == 1, 1, 'last'); % Find last bone element that correspond to the bone ST interface
end
%% 
profile3 = profile;
%% Plot map
figure;
for i = 1:snapNb
    snapName = sprintf('tx_%02d/T11_%03d.snp2D',tx, i);
    snap = SimSonic2DReadSnp2D(fullfile(simulationPath, snapName));
    
    % Plot
    imagesc(X, Z, snap.Data)
    hold on 
    plot(X, profile*grid.step, 'w', 'LineWidth', 2)
    hold on
    axis ij image, %shading interp
    pause(0.1) % Pause between the display of each snap
end

%%
figure
subplot(3,1,2)
p = plot(X, profile1*grid.step -10);
p.Color = "#0B5394";
p.LineWidth = 3;
axis equal
xlabel('Width (mm)', 'Interpreter', 'latex',  'FontSize', 22);
ylabel('Depth (mm)', 'Interpreter', 'latex',  'FontSize', 22);
ax = gca; 
ax.FontSize = 16;
legend(sprintf('%s = 0.38 mm, %s = 0.5 mm', texlabel('sigma'), texlabel('rho')) ...
    )
subplot(3,1,1)
p = plot(X, profile2*grid.step -10);
p.Color = "#0B5394";
p.LineWidth = 3;
axis equal
xlabel('Width (mm)', 'Interpreter', 'latex',  'FontSize', 22);
ylabel('Depth (mm)', 'Interpreter', 'latex',  'FontSize', 22);
legend(sprintf('%s = 0.38 mm, %s = 4.0 mm', texlabel('sigma'), texlabel('rho')))
ax = gca; 
ax.FontSize = 16;
title(sprintf('Profile regarding rms (%s) and correlation length(%s)', texlabel('rho'), texlabel('sigma')), 'FontSize', 24);
subplot(3,1,3)
p = plot(X, profile3*grid.step -10);
p.Color = "#0B5394";
p.LineWidth = 3;
xlabel('Width (mm)', 'Interpreter', 'latex',  'FontSize', 22);
ylabel('Depth (mm)', 'Interpreter', 'latex',  'FontSize', 22);
axis equal
legend(sprintf('%s = 0.08 mm, %s = 0.5 mm', texlabel('sigma'), texlabel('rho')))
ax = gca; 
ax.FontSize = 16;