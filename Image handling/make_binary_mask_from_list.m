function [path_to_mask] = make_binary_mask_from_list(list,filename,dims,dummy_image)
% Takes a voxel list and turns it into a binary image. This can then
% conveniently be compared to the subject's anatomy, for example with
% SPM's check reg function or mricro.
%
% Input arguments: 
% - list: an nx1 list of voxels (binary)
% - filanem: filepath (without file suffix) under which mask will be saved
% - dims: dimensions that the mask should have
% - dummy_image: dummy image that neads to be read in (e.g., a 'mask.img')
%
% - C. Utzerath 2014-15
%

%% Set input
fname = filename;
dims = dims; 
list = list;
example_mask = dummy_image;

%% Try reshaping
try
    reshaped_list = reshape(list,dims);
catch
    disp('Reshaping image went wrong - inconsistent image dimensions!')
end

%% Load in spm mask image and re-save with own data
% Read in dummy image
V = spm_vol(example_mask);

% Save your own data in the structure of the dummy image, under new fname
fname = [fname,'.nii'];
V.fname = fname;
spm_write_vol(V,reshaped_list);
disp(['Saved mask image to look at at: ' filename '.nii'])

path_to_mask = fname;

end

