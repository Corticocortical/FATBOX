function [freqs,amps] = make_frequency_spectrum(Y,FS,T,display)
% Creates a plottable frequency spectrum from an input signal.
%
%
% Note that this is not entirely the same as power density spectrum. Use
% make_power_density_spectrum for that instead. 
%
% Input arguments; 
% - Y: vector. The signal whose spectrum to plot
% - FS: scalar. Sampling frequency of the measured signal (Hertz)
% - T: vector. Time points of measured signal
% - display: binary. Plot the frequency spectrum?
%
% Output:
% - freqs: vector with different frequency components (Hertz), i.e. x-axis
% - amps: vector with amplitudes of corresponding frequency components
%
%% Create fake input arguments for manual mode
%https://dadorran.wordpress.com/2014/02/20/plotting-frequency-spectrum-using-matlab/

%% Do
fY = abs(ifft(Y));      % fft alternative A
N = numel(T);           % Length of signal

% Frequency axis
idx = 1:numel(fY) / 2;
frequencies = FS*(0:(N-1)) / N;

% Frequencies and amplitudes
freqs = frequencies(idx);
amps  = 2*fY(idx);

if display
    figure   
    plot(freqs,amps)
    xlabel('Frequency (cycles/second)'); ylabel('Amplitude');
    title('Frequency spectrum of signal')
end


end

