function [f,power] = make_power_density_spectrum(Y,FS,T, display)
% Creates a power density spectrum from a signal Y. Optionally plots
% spectrum.
%
% Requires MATLAB signal processing toolbox.
%
% 
% Input arguments; 
% - Y: vector. The signal whose spectrum to plot
% - FS: scalar. Sampling frequency of the measured signal (Hertz)
% - T: vector. Time points of measured signal
% - display: binary. Plot the frequency spectrum?
%
% Output arguments:
% - power: vector. Power of each sampled frequecy.
% - f: vector. Sampled frequencies.
%
%
%% Do

[power,f] = periodogram(Y,[],length(Y),FS);

if display
    figure
    plot(f,power)
    xlabel('Frequencies','fontsize',12)
    ylabel('Power (spectral density)','fontsize',12)
    title('Power spectrum of signal','fontsize',12)
    box off
end

end

