function FullImage = stitching_TBM(ImageTissue,ImageBone,ImageMarrow,P_Periosteum,P_Endosteum,X,Z,wavelength_skin,wavelength_bone)

    pixel_size = Z(2)-Z(1);

    depth_periosteal_surf = P_Periosteum(1)*X.^2+P_Periosteum(2)*X+P_Periosteum(3);
    for ll=1:length(X)
        [Y,Ind] = min(abs(Z-depth_periosteal_surf(ll)));
        new_index_depth_periosteal_surf(ll) = Ind;
        index_stitching = new_index_depth_periosteal_surf + round(wavelength_skin/pixel_size);
    end

    depth_endosteal_surf = P_Endosteum(1)*X.^2+P_Endosteum(2)*X+P_Endosteum(3);
    for ll=1:length(X)
        [Y,Ind] = min(abs(Z-depth_endosteal_surf(ll)));
        new_index_depth_endosteal_surf(ll) = Ind;
        index_stitching_endosteum = new_index_depth_endosteal_surf + round(wavelength_skin/pixel_size);
    end

    nb_pts_smooth_stitching = round(wavelength_bone/pixel_size);
    nb_pts_smooth_stitching = uint32(nb_pts_smooth_stitching + (nb_pts_smooth_stitching/2~=round(nb_pts_smooth_stitching/2)));
    win = tukeywin(2*nb_pts_smooth_stitching+1,1);
    
    win_fade_down = win(nb_pts_smooth_stitching+1:2*nb_pts_smooth_stitching+1);
    win_fade_down = repmat(win_fade_down,1,length(X));
    win_fade_up = win(1:nb_pts_smooth_stitching+1);
    win_fade_up = repmat(win_fade_up,1,length(X));
    
    
    mask_stitching_I1 = ones(length(Z),length(X));
    for ii=1:length(X)
        mask_stitching_I1(index_stitching(ii)-nb_pts_smooth_stitching/2:index_stitching(ii)+nb_pts_smooth_stitching/2,ii) = win_fade_down(:,ii);
        mask_stitching_I1(index_stitching(ii)+nb_pts_smooth_stitching/2+1:end,ii) = 0;
    end
    
    mask_stitching_I2 = ones(length(Z),length(X));
    for ii=1:length(X)
        mask_stitching_I2(1:index_stitching(ii)-nb_pts_smooth_stitching/2-1,ii) = 0;
        mask_stitching_I2(index_stitching(ii)-nb_pts_smooth_stitching/2:index_stitching(ii)+nb_pts_smooth_stitching/2,ii) = win_fade_up(:,ii);
        mask_stitching_I2(index_stitching_endosteum(ii)-nb_pts_smooth_stitching/2:index_stitching_endosteum(ii)+nb_pts_smooth_stitching/2,ii) = win_fade_down(:,ii);
        mask_stitching_I2(index_stitching_endosteum(ii)+nb_pts_smooth_stitching/2+1:end,ii) = 0;
    end
    
    mask_stitching_I3 = ones(length(Z),length(X));
    for ii=1:length(X)
        mask_stitching_I3(1:index_stitching_endosteum(ii)-nb_pts_smooth_stitching/2-1,ii) = 0;
        mask_stitching_I3(index_stitching_endosteum(ii)-nb_pts_smooth_stitching/2:index_stitching_endosteum(ii)+nb_pts_smooth_stitching/2,ii) = win_fade_up(:,ii);
    end

    mask_stitching_I1 = mask_stitching_I1(1:length(Z),:);
    mask_stitching_I2 = mask_stitching_I2(1:length(Z),:);
    mask_stitching_I3 = mask_stitching_I3(1:length(Z),:);

    
    
    tmp1 = ImageTissue.*mask_stitching_I1;
    tmp2 = ImageBone.*mask_stitching_I2;
    tmp3 = ImageMarrow.*mask_stitching_I3;
    FullImage = tmp1 + tmp2 + tmp3;

   
    return