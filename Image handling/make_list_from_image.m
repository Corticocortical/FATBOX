function [list,dims] = make_list_from_image(imagepath)
% Reshapes a NIFTI-image into a nvox X 1 list. Useful for binary mask
% images or to sort images easily.
%
% Uses SPM function to read images.
%
% C. Utzerath, 2014-15

%% Read image
volume = spm_vol(imagepath);
image  = spm_read_vols(volume);

%% Reshape image, return
dims = volume.dim;
list = reshape(image,max(cumprod(dims)),1);

end

