function [ray_specu] = get_tof_specu_two_layers(P0,T,interface_param,specular_object,V1,V2)
    % T et R son confondus
    %     xr = T(1);zr = T(2);
    xp0 = P0(1);zp0 = P0(2);
    xt0 = T(1);zt0 = T(2);    
    a_p = specular_object(1);b_p = specular_object(2);
    c_p = zp0-a_p*xp0^2-b_p*xp0;
    a_i = interface_param(1);b_i = interface_param(2);c_i = interface_param(3);
    lower_bound = -5e-2;upper_bound = -lower_bound;
    alpha = @ (xt,zt,xp,zp) atan((xp-xt)./(zp-zt));
    thetaL = @(xi) -atan(2*a_i*xi+b_i);
    z_i = @(xi) a_i*xi^2+b_i*xi+c_i;
    zp = @(xp) a_p*xp^2+b_p*xp+c_p;
    x_i = @(xp) get_incidence_point(xp,zp(xp),interface_param, specular_object,lower_bound, upper_bound);
    f_specular_refra =@ (xp,xi) (1/V1*sin(alpha(xt0,zt0,xi,z_i(xi))-thetaL(xi))...
        -1/V2*sin(alpha(xi,z_i(xi),xp,zp(xp))-thetaL(xi)));
    f_specular_refra_fast_abs = @ (xp) abs(f_specular_refra(xp,x_i(xp)));
    f_specular_refra_fast = @ (xp) (f_specular_refra(xp,x_i(xp)));

%     rng default % For reproducibility
%     opts = optimoptions(@fmincon,'Algorithm','sqp');
%     problem = createOptimProblem('fmincon','objective',...
%         f_specular_refra_fast,'x0',0,'lb',lower_bound,'ub',upper_bound,'options',opts);
%     gs = GlobalSearch;

%     x_specu = fminbnd(f_specular_refra_fast_abs,lower_bound, upper_bound);
    x_specu= fzero(f_specular_refra_fast,0e-3);
%     [x_specu,f] = run(gs,problem);    
    z_specu = zp(x_specu);
    xi_specu = x_i(x_specu);
    zi_specu = z_i(x_specu);
    ray1 = 2*hypot(xt0-xi_specu,zt0-zi_specu);
    ray2 = 2*hypot(x_specu-xi_specu,z_specu-zi_specu);
    ray_specu = ray1/V1 + ray2/V2;
%     figure, fplot(f_specular_refra_fast,[lower_bound upper_bound])
%     hold on
%     plot(x_specu,f_specular_refra_fast(x_specu),'r.')
% %     plot(x_specu_0,f_specular_refra_fast(x_specu_0),'b.')
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