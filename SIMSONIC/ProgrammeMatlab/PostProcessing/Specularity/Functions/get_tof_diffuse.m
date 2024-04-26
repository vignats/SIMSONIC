function [ray_specu] = get_tof_diffuse(P0,T,periosteum,V1,V2)
    % No lens
    % T et R son confondus
    %     xr = T(1);zr = T(2);
    xp0 = P0(1);zp0 = P0(2);
    xt0 = T(1);zt0 = T(2);    
    a_i = periosteum(1);b_i = periosteum(2);c_i = periosteum(3);
    lower_bound = -5e-2;upper_bound = -lower_bound;
    alpha = @ (xt,zt,xp,zp) atan((xp-xt)./(zp-zt));
    thetaL = @(xi) -atan(2*a_i*xi+b_i);
    z_i = @(xi) a_i*xi^2+b_i*xi+c_i;
    f_specular_refra =@ (xi) (1/V1*sin(alpha(xt0,zt0,xi,z_i(xi))-thetaL(xi))...
        -1/V2*sin(alpha(xi,z_i(xi),xp0,zp0)-thetaL(xi)));
    % f_specular_refra_fast_abs = @ (xp) abs(f_specular_refra(xp,x_i(xp)));
    % f_specular_refra_fast = @ (xp) (f_specular_refra(xp,x_i(xp)));
    xi= fzero(f_specular_refra,0e-3);
    zi= z_i(xi);
    ray1 = 2*hypot(xt0-xi,zt0-zi);
    ray2 = 2*hypot(xp0-xi,zp0-zi);
    ray_specu = ray1/V1 + ray2/V2;
end
