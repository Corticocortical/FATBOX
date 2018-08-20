function [tvec] = read_spm_tmap(fpath)
% Reads an SPM tmap and returns it in vectorized format along with its
% dimensions.
%
% Input argument:
%   - fpath: filename of tmap
% 
% Output: 
%   - tvec: Nx1 vector of the map
%   - dims: dimensions of the original map
%
% C. Utzerath 2014-15

%% Read map
vol  = spm_vol(fpath);
img  = spm_read_vols(vol);
dims = vol.dim;
tvec = reshape(img,numel(img),1);

end

