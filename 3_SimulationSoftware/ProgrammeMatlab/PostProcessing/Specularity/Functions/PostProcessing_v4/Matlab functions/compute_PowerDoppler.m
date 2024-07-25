function [PDim,PDim_SpatialAv]=compute_PowerDoppler(IQ_im,Average_kernel_size_x,Average_kernel_size_z,kernel_size_temporal_smoothing)


%tmp = squeeze(IQ_im.*conj(IQ_im));
tmp = squeeze(abs(IQ_im).^2);

if (Average_kernel_size_z>1 && Average_kernel_size_x>1)
    %-----------------------------------------%
    %---------- Spatial smoothing ------------%
    %-----------------------------------------%
    sAvg = ones(Average_kernel_size_z,Average_kernel_size_x)./Average_kernel_size_x/Average_kernel_size_z;
    for it=1:size(tmp,3)    
        tmp(:,:,it) = filter2(sAvg,tmp(:,:,it)); 
    end
end


if (kernel_size_temporal_smoothing>1)
    %-----------------------------------------%
    %--------- Temporal smoothing ------------%
    %-----------------------------------------%
    ht=ones(1,kernel_size_temporal_smoothing)/kernel_size_temporal_smoothing;
    tmp = filter(ht,1,tmp,[],3);
end

PDim = tmp;
PDim_SpatialAv = squeeze(mean(PDim,[1,2]));

