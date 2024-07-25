function [stitching_mask_Tissue, stitching_mask_Bone]=get_stitching_masks...
    (X,Z, peribola_coef, center_freq,V1,V2)

    pixel_size = Z(2)-Z(1);
    lambda_soft = V1/center_freq; lambda_bone = V2/center_freq;
    Nz = numel(Z); Nx =numel(X);
    depth_start_periosteal_surf = peribola_coef(1)*X.^2+peribola_coef(2)*X+peribola_coef(3);
    
    for ll=1:Nx
        [~,Ind] = min(abs(Z-depth_start_periosteal_surf(ll)));
        new_index_depth_start_periosteal_surf(ll) = Ind;
        index_stitching = new_index_depth_start_periosteal_surf + round(lambda_soft/pixel_size);
    end
    
    nb_pts_smooth_stitching = round(lambda_bone/pixel_size);
    nb_pts_smooth_stitching = uint32(nb_pts_smooth_stitching + (nb_pts_smooth_stitching/2~=round(nb_pts_smooth_stitching/2)));
    win = tukeywin(2*nb_pts_smooth_stitching+1,1);
    
    win_I0 = repmat(win(nb_pts_smooth_stitching+1:2*nb_pts_smooth_stitching+1),1,Nx);
    stitching_mask_Tissue = ones(Nz,Nx);
    for ii=1:Nx
        stitching_mask_Tissue(index_stitching(ii)-nb_pts_smooth_stitching/2:index_stitching(ii)+nb_pts_smooth_stitching/2,ii) = win_I0(:,ii);
        stitching_mask_Tissue(index_stitching(ii)+nb_pts_smooth_stitching/2+1:end,ii) = 0;
    end
    
    
    win_I2 = repmat(win(1:nb_pts_smooth_stitching+1),1,Nx);
    stitching_mask_Bone = ones(Nz,Nx);
    for ii=1:Nx
        stitching_mask_Bone(index_stitching(ii)-nb_pts_smooth_stitching/2:index_stitching(ii)+nb_pts_smooth_stitching/2,ii) = win_I2(:,ii);
        stitching_mask_Bone(1:index_stitching(ii)-nb_pts_smooth_stitching/2-1,ii) = 0;
    end
    
    stitching_mask_Tissue = stitching_mask_Tissue(1:Nz,:);
    stitching_mask_Bone = stitching_mask_Bone(1:Nz,:);
    % figure, imagesc(mask_stitching_Tissue)
    % figure, imagesc(mask_stitching_Bone)
end