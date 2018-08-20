function [] = reslice_image(image_to_reslice,image_defining_space)
% Reslices an image so it has the same dimensions as the other. Images
% should already be warped into the same template and should be properly
% coregistered. Uses 4th B-Spline interpolation.
% 
% Places the resliced image in the same folder as the image to reslice,
% prefixed with an 'r'.
%
% Input argument:
% - image_to_reslice: string. Path to imgae.
% - image_defining_space: string. Path to the image that you reslice
%   to.
%
% C. Utzerath, 2015

%% Defne batch
matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.write.ref = {image_defining_space};
matlabbatch{1}.spm.spatial.coreg.write.source = {image_to_reslice};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

%% Run
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
disp('Resliced.')

end

