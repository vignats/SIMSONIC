function [ray_specu] = get_tof_specu_one_layers(P0,T,specular_object,V0)
    % T et R son confondus
    %     xr = T(1);zr = T(2);
    xp0 = P0(1);zp0 = P0(2);
    xt0 = T(1);zt0 = T(2);    
    a_p = specular_object(1);b_p = specular_object(2);
    c_p = zp0-a_p*xp0^2-b_p*xp0;
    lower_bound = -5e-2;upper_bound = -lower_bound;

    zp = @(xp) a_p*xp^2+b_p*xp+c_p;

    f_specu_simpl = @ (xp) atan((xt0-xp)+(zt0-zp(xp))*(2*a_p*xp+b_p));
    x_specu = fzero(f_specu_simpl,0e-3);

    z_specu = zp(x_specu);
    ray_specu = 2*hypot(xt0-x_specu,zt0-z_specu)/V0;

    show=0;
    if show
        figure, 
        X = lower_bound:1e-3:upper_bound;
        plot(X*1e3,polyval([a_p b_p c_p],X)*1e3),
        axis ij,
        ylim([0 20])
        hold on
    %     plot(X*1e3,polyval([a_p b_p c_p],X)*1e3),
        plot(X*1e3,polyval([a_i b_i c_i],X)*1e3),
        plot(xp0*1e3,zp0*1e3,'kx')
        plot(xt0*1e3,0*1e3,'o')
    %     plot(x_specu*1e3,z_specu*1e3,'+')
        plot(xi_specu*1e3,zi_specu*1e3,'v')
        plot([xi_specu xt0]*1e3,[zi_specu zt0]*1e3,'--')
        plot([xi_specu x_specu]*1e3,[zi_specu z_specu]*1e3)
    end
end
function [xi] = get_incidence_point(xp,zp,interface_param, specular_object,lower_bound,upper_bound)
    a_p = specular_object(1);b_p = specular_object(2);
%     c_p = zp0-a_p*xp^2-b_p*xp;
%     a_i = interface_param(1);b_i = interface_param(2);c_i = interface_param(3);    
%     lower_bound = -5e-2;upper_bound = -lower_bound;
    z_i = @(x) polyval(interface_param,x);
    f_specu_abs = @ (xi) abs(atan((xi-xp)/(z_i(xi)-zp))+atan(2*a_p*xp+b_p));
    f_specu = @ (xi) atan((xi-xp)/(z_i(xi)-zp))+atan(2*a_p*xp+b_p);
    f_specu_simpl = @ (xi) atan((xi-xp)+(z_i(xi)-zp)*(2*a_p*xp+b_p));
%         xi=fzero(f_specu,0e-3);
        xi = fzero(f_specu_simpl,0e-3);
%     figure, 
%     subplot 211 
%     fplot(f_specu_abs,[lower_bound upper_bound]) 
%     hold on
%     fplot(f_specu_simpl,[lower_bound upper_bound],'k--')
%     subplot 212
%     fplot(f_specu,[lower_bound upper_bound])
%     hold on
% 
% % %     plot(x0,f_specu(x0),'r.')
%     plot(xi,f_specu(xi),'kx')
%     f_specu_abs(xi)

end