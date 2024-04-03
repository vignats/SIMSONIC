function SPECU_MODEL_MATLAB_BONE=function_get_pix_simplified_tissue_model(xp,zp,interface_parameters,SPEED,XR,...
    excitation_signal,FS)
    NTX = numel(XR);
    if numel (interface_parameters)==1
        a = -tand(interface_parameters); b = zp-a*xp;
        SP_PARABOLA = [0 a b];
    elseif numel(interface_parameters)==2
        c = zp-interface_parameters(1)*xp^2-interface_parameters(2)*xp;
        SP_PARABOLA = [interface_parameters(1) interface_parameters(2) c];
    end
        P0 = [xp,zp];
        TOF_MATLAB_BONE = zeros([1 NTX]);
        TOF_DIFF_BONE = zeros([1 NTX]);
        for tx=1:NTX
            T = [XR(tx) 0];
            TOF_MATLAB_BONE(tx)=get_specular_time_of_flight(P0,T,SP_PARABOLA,SPEED);
            TOF_DIFF_BONE(tx) = 2*hypot(xp-T(1),zp)/SPEED;
        end
        DELTA_T = round((TOF_DIFF_BONE-TOF_MATLAB_BONE+0.6e-6)*FS);
        DELTA_T(DELTA_T<=0) = nan;
        SPECU_MODEL_MATLAB_BONE = zeros(size(TOF_MATLAB_BONE));
        SPECU_MODEL_MATLAB_BONE(DELTA_T<=numel(excitation_signal)) = excitation_signal(DELTA_T(DELTA_T<=numel(excitation_signal)));
        SPECU_MODEL_MATLAB_BONE = SPECU_MODEL_MATLAB_BONE/max(SPECU_MODEL_MATLAB_BONE);
end

function [ray_specu] = get_specular_time_of_flight(P0,T,specular_object,SPEED)
    % T et R son confondus
    %     xr = T(1);zr = T(2);
    xp0 = P0(1);zp0 = P0(2);
    xt0 = T(1);zt0 = T(2);    
    a_p = specular_object(1);b_p = specular_object(2);
    c_p = zp0-a_p*xp0^2-b_p*xp0;
    zp = @(xp) a_p*xp^2+b_p*xp+c_p;
    f_specu_simpl = @ (xp) atan((xt0-xp)+(zt0-zp(xp))*(2*a_p*xp+b_p));
    x_specu = fzero(f_specu_simpl,0e-3);

    z_specu = zp(x_specu);
    ray_specu = 2*hypot(xt0-x_specu,zt0-z_specu)/SPEED;


end