function [mean_vertical_spatial_period_image,median_vertical_spatial_period_image,mean_lateral_spatial_period_image,median_lateral_spatial_period_image,spatial_period]=estimate_spatial_period_image(IQim,X,Z,assumed_spatial_period_image,mask_ROI,diplay_YES)

mask_NaN = repmat(mask_ROI,1,1,size(IQim,3));

%------------------------------------------------%
%------- estimation spatial period in image -----%
%----------------- ALL FRAMES -------------------%
%------------------------------------------------%

% Loupas's method estimates the vertical spatial period in image
pixel_size = Z(2) - Z(1);
tmp = angle(IQim(2:end,:,:).*conj(IQim(1:end-1,:,:)));
exact_spatial_period_image = 1./(tmp/2/pi)*pixel_size;
exact_spatial_period_image = abs(exact_spatial_period_image);
exact_spatial_period_image(~mask_NaN(1:end-1,:,:)) = NaN;
exact_spatial_period_image(exact_spatial_period_image>assumed_spatial_period_image*5) = NaN; 
mean_vertical_spatial_period_image = mean(exact_spatial_period_image(:),'omitnan');
median_vertical_spatial_period_image = median(exact_spatial_period_image(:),'omitnan');

if diplay_YES
figure(123)
subplot 121
histogram(exact_spatial_period_image(~isnan(exact_spatial_period_image))*1e3,100)
title(['vertical spatial period in image [mm]'])
set(gca,'fontsize',14)

figure(22)
subplot 131
imagesc(X*1e3,Z*1e3,real(IQim(:,:,end/2)))
ylabel('depth (mm)')
xlabel('width (mm)')
axis image
set(gca,'fontsize',14)
subplot 132
imagesc(X*1e3,Z(1:end-1)*1e3,median(exact_spatial_period_image*1e3,3,'omitnan'))
ylabel('depth (mm)')
xlabel('width (mm)')
colorbar
axis image
title(['vertical spatial period in image [mm]'])
set(gca,'fontsize',14)
end

% Loupas's method estimates the lateral spatial period in image
pixel_size = X(2) - X(1);
tmp = angle(IQim(:,2:end,:).*conj(IQim(:,1:end-1,:)));
exact_spatial_period_image = 1./(tmp/2/pi)*pixel_size;
exact_spatial_period_image = abs(exact_spatial_period_image);
exact_spatial_period_image(~mask_NaN(:,1:end-1,:)) = NaN;
exact_spatial_period_image(exact_spatial_period_image>assumed_spatial_period_image*20) = NaN; 
mean_lateral_spatial_period_image = mean(exact_spatial_period_image(:),'omitnan');
median_lateral_spatial_period_image = median(exact_spatial_period_image(:),'omitnan');

spatial_period = 1/sqrt(median_vertical_spatial_period_image^(-2) + median_lateral_spatial_period_image^(-2));

if diplay_YES
figure(123)
subplot 122
histogram(exact_spatial_period_image(~isnan(exact_spatial_period_image))*1e3,100)
title(['lateral spatial period in image [mm]'])
set(gca,'fontsize',14)

figure(22)
subplot 133
imagesc(X*1e3,Z(1:end-1)*1e3,median(exact_spatial_period_image*1e3,3,'omitnan'))
ylabel('depth (mm)')
xlabel('width (mm)')
colorbar
axis image
title(['lateral spatial period in image [mm]'])
set(gca,'fontsize',14)
end

end