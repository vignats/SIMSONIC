function [recorded] = LoadRfData(probe, simu_dir)
% LoadRfData generates matrix of the RF data from all receiver transducer.
%
% Syntax:
%   rf = LoadRfData(probe, simu_dir)
%
% Input Arguments:
%   - probe: Structure containing information about the ultrasound probe.
%   - simu_dir : Path to the directory that contain the simulation
%   information for each transducer
%
% Output Arguments:
%   - rf : matrix of dimension time x receiver x transmitter
    for tx = 1 : probe.Nelements
        rf_path = fullfile(simu_dir, sprintf('tx_%02d/T11_main_T11.rcv2D',tx));
        signals = SimSonic2DReadRcv2D(rf_path);
        rf(:,:,tx) = signals.Signals;
    end 

    recorded = signals;
    recorded.Signals = rf;
end 
