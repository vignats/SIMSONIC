function SPECU_MODEL_MATLAB_BONE=get_pix_simplified_bone_model_NEW(ix,iz,interface_parameters,geometry,XR,...
    excitation_signal,transmit_time,FS)
%     transmit_time = TOF.Time_T(:,iz,ix);%rec_time = TOF.Time_R(:,iz,ix);
%     transmit_angle_px=AOV_px.Angle_T(:,iz,ix);rec_angle_px=AOV_px.Angle_R(:,iz,ix);
%     transmit_angle_elem=AOV_elem.Angle_T(:,iz,ix);rec_angle_elem=AOV_elem.Angle_R(:,iz,ix);
    zp=geometry.Z(iz);xp = geometry.X(ix);
    lt = geometry.lens_thick;
    NTX = numel(XR);
    C_LENS = geometry.V0; C_TISSUE = geometry.V1; C_BONE= geometry.V2;
    C= [C_LENS C_TISSUE C_BONE 0 0];
        % [c_lens c_tissue c_bone aniso_coef aniso_shape_coef];
    if numel (interface_parameters)==1
        a = -tand(interface_parameters); b = zp-a*xp;
        SP_PARABOLA = [0 a b];
    elseif numel(interface_parameters)==2
        c = zp-interface_parameters(1)*xp^2-interface_parameters(2)*xp;
        SP_PARABOLA = [interface_parameters(1) interface_parameters(2) c];
    end
    PARABOLA_COEFS_PERI = geometry.peribola_coef;
    P = [PARABOLA_COEFS_PERI SP_PARABOLA lt];% for iz=1:NZ
%     SPECU_MODEL_MATLAB_BONE = zeros([1 NTX]);
%     if all(transmit_time)
        P0 = [xp,zp];
        TOF_MATLAB_BONE = zeros([1 NTX]);
        TOF_DIFF_BONE = zeros([1 NTX]);
        for tx=1:NTX
            T = [XR(tx) 0];
            TOF_MATLAB_BONE(tx)=get_tof_specu_two_layers(P0,T,P(1:3),P(4:6),C(2),C(3));
            TOF_DIFF_BONE(tx) = get_tof_diffuse(P0,T,P(1:3),C(2),C(3));%2*(transmit_time(tx));
        end
        DELTA_T = round((TOF_DIFF_BONE-TOF_MATLAB_BONE+0.6e-6)*FS);
        DELTA_T(DELTA_T<=0) = nan;
        SPECU_MODEL_MATLAB_BONE = zeros(size(TOF_MATLAB_BONE));
        SPECU_MODEL_MATLAB_BONE(DELTA_T<=numel(excitation_signal)) = excitation_signal(DELTA_T(DELTA_T<=numel(excitation_signal)));
        SPECU_MODEL_MATLAB_BONE = SPECU_MODEL_MATLAB_BONE/max(SPECU_MODEL_MATLAB_BONE);
%     end
% TOF_DIFF_BONE - transmit_time.'
end