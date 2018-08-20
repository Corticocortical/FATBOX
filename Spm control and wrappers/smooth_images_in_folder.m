function [ output_args ] = smooth_images_in_folder(folder,prefix,kernel)
% Smooths functional images of type <filter> that are foud in <folder>,
% with a smoothing kernel of <kernel> mm.
%
% Input arguments:
% - folder: string. folder in which the files are located.
% - prefix: string. file prefix to scan for, e.g., '' 
% - kernel: 1x3 vec with resolution, e.g., [8 8 8].
%
% Created for SPM12.
% % C. Utzerath, 2014-16
%
%% Parse example input
folder = folder;
prefix = prefix;
kernel = kernel;

%% Make batch
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {folder};
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = prefix;
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('File Selector (Batch Mode): Selected Files (arf*.nii)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = kernel;
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';

%% Run batch
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
disp(['Smoothed images at ' folder])

end

