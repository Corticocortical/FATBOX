function [freq] = convert_rad_per_sample_to_frequency(rad_per_sample,SF)
% Converts radians per sample, as MATLAB filters use them, to frequencies.
%
% Input arguments:
% - rad_per_sample: scalar/vector. The values to be converted. If supplied
%                   as vector, script will convert each element.
% - SF. Scalar. Sampling frequency.
%
% Output arguments:
% - freq: sclar/vector. Frequencies corresponding to input. 
%
%
%% Do
T = 1/SF; % Sampling rate to sampling interval

freq = [];
for i = 1:numel(rad_per_sample)
    freq(i) =  rad_per_sample(i) / (2*pi*T);
end

end

