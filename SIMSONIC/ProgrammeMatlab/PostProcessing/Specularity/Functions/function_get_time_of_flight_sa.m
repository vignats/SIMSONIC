function [Time_T, Time_R] = function_get_time_of_flight_sa(X,Z,sos,XT,ZT,XR,ZR)
    NTx = numel(XT);NRx = numel(XR);
    NZ = numel(Z);NX = numel(X);
    Time_T = zeros([NTx NZ NX]);
    Time_R = zeros([NRx NZ NX]);
    iz=0;
    for z=Z
        iz=iz+1;
        ix=0;
        for x=X
            ix=ix+1;
            tx_flight = hypot(x-XT,z-ZT)./sos;
            rx_flight = hypot(x-XR,z-ZR)./sos;
            Time_T(:,iz,ix) = tx_flight;
            Time_R(:,iz,ix) = rx_flight;
        end
    end
end