function [Angle_T,Angle_R] = function_get_angle_of_view_sa(X,Z,XT,ZT,XR,ZR)
    NTx = numel(XT);NRx = numel(XR);
    NZ = numel(Z);NX = numel(X);
    Angle_T  = zeros([NTx NZ NX]);Angle_R = zeros([NRx NZ NX]);
    zi = 0;
    for z=Z
        zi=zi+1;xi=0;        
        for x=X
            xi=xi+1;
            Angle_T(:,zi,xi) = atand((x-XT)./(z-ZT));
            Angle_R(:,zi,xi) = atand((x-XR)./(z-ZR));
        end
    end

end