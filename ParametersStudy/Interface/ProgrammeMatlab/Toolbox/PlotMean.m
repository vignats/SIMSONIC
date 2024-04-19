linearROI_RMS_05 = zeros(numel(rmsAll), length(reconstruction.Xmm));
for i = 1:numel(rmsAll)
    linearROI_RMS_05(i, :) = SpecularAll.SpecuProba{i,1}{1}.linearROI;
    linearROI_RMS_1(i, :) = SpecularAll.SpecuProba{i,2}{1}.linearROI;
    linearROI_RMS_2(i, :) = SpecularAll.SpecuProba{i,3}{1}.linearROI;
    linearROI_RMS_4(i, :) = SpecularAll.SpecuProba{i,4}{1}.linearROI;
end

figure;
subplot(2,2,1)
surf(reconstruction.Xmm, rmsAll, linearROI_RMS_05, 'EdgeColor', 'none');
xlabel('Lateral Position (mm)');
ylabel('Root mean square (mm)');
zlabel('Mean Specular Probability ');
title('Mean Specular Probability for a correlation length of 0.5mm');
colorbar

subplot(2,2,2)
surf(reconstruction.Xmm, rmsAll, linearROI_RMS_1, 'EdgeColor', 'none');
xlabel('Lateral Position (mm)');
ylabel('Root mean square (mm)');
zlabel('Mean Specular Probability ');
title('Mean Specular Probability','Correlation length of 1mm');
colorbar

subplot(2,2,3)
surf(reconstruction.Xmm, rmsAll, linearROI_RMS_2, 'EdgeColor', 'none');
xlabel('Lateral Position (mm)');
ylabel('Root mean square (mm)');
zlabel('Mean Specular Probability ');
title('Mean Specular Probability','Correlation length of 2mm');
colorbar

subplot(2,2,4)
surf(reconstruction.Xmm, rmsAll, linearROI_RMS_4, 'EdgeColor', 'none');
xlabel('Lateral Position (mm)');
ylabel('Root mean square (mm)');
zlabel('Mean Specular Probability ');
title('Mean Specular Probability','Correlation length of 4mm');
colorbar
