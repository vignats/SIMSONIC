function rx_idx=get_idx_of_closest_value_v2(angle_r_deg,specular_angle)
    NRx = size(angle_r_deg,1);
    c = angle_r_deg(end:-1:1);
    bin      = [-Inf; (c(1:end-1) + c(2:end)) * 0.5; Inf];
    try
        rx_idx = discretize(specular_angle, bin);
        rx_idx = NRx+1-rx_idx;
    catch err
        if strcmp(err.identifier,'MATLAB:discretize:InvalidSecondInput')
            [cs,idx_sorted] = sort(c);
            bin      = [-Inf; (cs(1:end-1) + cs(2:end)) * 0.5; Inf];
            rx_idx = discretize(specular_angle, bin);
            rx_idx = idx_sorted(rx_idx)';
        else
            err.rethrow;
        end
    end
end