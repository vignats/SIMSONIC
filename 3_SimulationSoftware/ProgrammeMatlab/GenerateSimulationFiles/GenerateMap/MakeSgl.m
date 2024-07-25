function[] = MakeSgl(param, grid, medium, signal, print, simu_dir) 
% MakeSgl generates the excitation signal and writes it to Signal.sgl file in the simulation directory.
%
% Syntax:
%   MakeSgl(param, grid, medium, signal, simu_dir, print)
%
% Input Arguments:
%   - param: Structure containing simulation parameters.
%   - grid: Structure containing information about the grid dimensions.
%   - medium: Structure containing information about the medium properties.
%   - signal: Structure containing information about the excitation signal.
%   - print: Logical value indicating whether to plot the excitation signal.
%   - simu_dir: Directory where the simulation files are saved, 
%           if not indicated, the signal is only computed on matlab, 
%           not registered in the simulation directory.
 

    std_t = 1.52/(2*pi*signal.B);                           % Temporal standard deviation, exact
    periode = signal.nb_periode/signal.fc;    
    
    % Time vector
    dt = param.cfl*grid.step*1e-3/(sqrt(2)*max(medium.cp));      % Time step (s)
    Fs = 1/dt;                                              % Sampling frequency (Hz)
    N = round(Fs * periode);                                % Number of sample
    t = linspace(0, periode, N);                            % Time vector corresponding to the frequency domain
    
    % Computation of the signal
    signal_out = exp(-(t-periode/2).^2/(2*std_t.^2)).*sin(2*pi*signal.fc*t);
    
    % Plot the signal in the time domain if print is true
    % Plot the signal in the time domain if print is true
    if print
        figure;
        plot(t, signal_out);
        xlabel('Time (us)');
        ylabel('Amplitude');
        title('Excitation signal in time domain');
    end

    % Execute the last line only if simu_dir is specified
    if nargin ==6 
        fprintf(['\n--- Signal saved in ', simu_dir(46:end-1)]);
        SimSonic2DWriteSgl(transpose(signal_out), fullfile(simu_dir, 'Signal.sgl'));
    end
end