%% Visualize signal patwhay 
X = (0:2999)*10e-6; X = X - mean(X); % X corresponds to size_2 and Z to size_1
Z = (0:1499)*10e-6;
Xmm = X*1e3;Zmm = Z*1e3;
figure
for k=1:40
    snap_name = sprintf('/home/salomevignat/Documents/SIMSONIC/Simulation/simulation_rough_interface/tx_01/T11_%03d.snp2D',k);
    snap = SimSonic2DReadSnp2D(snap_name);
    pcolor(Xmm, Zmm, snap.Data)
    axis ij image, shading interp
    pause(0.1)
end
