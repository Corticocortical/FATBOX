function [amplitudes] = qc_amplitudes_slices(volume)
% Checks a a volume's slices' deviation from the mean. Useful for detecting
% overactive slices.
% 
% Initially this script was plotting the slice mean, but Floris mentioned
% that it would be more meaningful to instead plot the deviation from the
% mean. This is done simply by subtracting each slice's mean from the slice
% voxels prior to the resulting matrix of deviations.
%
%
% Input argument: 
% - volume: string. Path to volume. 
%
% Output:
% - Amplitudes: vector containing mean amplitude per slice. 
%
% NEEDS: 
% SPM to read volumes. 
%
%
% Christian Utzerath 2014-15 (Donders Institute)


%% Do 
% Read slice
volume    = spm_vol(volume);
volume    = spm_read_vols(volume);

% Per slice, calculate each voxels' deviation from grand mean
newvolume = [];
for slice = 1:size(volume,3)
    slicemat             = volume(:,:,slice);
    slicemean            = mean(mean(slicemat));
    newvolume(:,:,slice) = slicemat - slicemean;
end

% Compute average
av = squeeze(mean(newvolume,1)); % average away first dimension
av = squeeze(mean(av,1));     % average away remaining dimension


% Clear, output
clear volume
amplitudes = av;
clear av;

end

