function [rad_per_sample] = convert_freq_to_rad_per_sample(frequency,SF)
% Converts a frequency to radians per sample, which is used by MATLAB's
% filtering tools. 
%
% Input argument: 
% - frequency. Scalar. The frequency to convert. If frequency is supplied
%              as vector, conversion will be performed for each element
%              of it.
% - SF. Scalar. Sampling rate (events/second).
%
% Output arguments:
% - rad_per_sample: scalar. Radians per sample. If frequency is supplied
%                   as vector, conversion will be performed for each element
%                   of it. 
%
%% Do
T = 1/SF;                                       % Compute sampling interval from sampling rate

rad_per_sample = [];
for i = 1:numel(frequency)
    rad_per_sample(i) =  2 *  pi * frequency(i) * T;      % Cf. sinus formulas: 2*pi*f*t = f
end




end

