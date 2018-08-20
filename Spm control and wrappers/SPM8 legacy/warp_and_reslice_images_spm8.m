function[] =  warp_and_reslice_images_spm8(params,vx_dim,images,spaceImage)
% Uses SPM to bring a set of images into another space by warping them into
% that space and then reslicing.
%
% All involved images should already be properly co-registered.
% Normalization parameters can best be obtained by segmenting a subject
% (does not take long).
%
% Input arguments:
% - params: string. file path to param file.
% - vx_dim: 1x3 vector. Voxel dimensions of output space
% - spaceImage: string. file path to image defining space.
% - images: cell with strings that contain file paths to images to process.

% params = '/home/predatt/chrutz/ImpexpII/Data/S14/fMRI/Anatomical/s150504084444DST131221107523026393-0002-00001-000176-01_seg_inv_sn.mat';
% spaceImage = '/home/predatt/chrutz/ImpexpII/Data/S14/fMRI/Results/Pilot pipeline 1/RSA pipeline 1/RSA/localizer GLM/mask.img';
% images     =     {
%                 '/home/predatt/chrutz/ImpexpII/Data/Group fMRI/Group masks/rFG1_l.nii'
%                 '/home/predatt/chrutz/ImpexpII/Data/Group fMRI/Group masks/rFG1_r.nii'
%                 '/home/predatt/chrutz/ImpexpII/Data/Group fMRI/Group masks/rFG2_l.nii'
%                 '/home/predatt/chrutz/ImpexpII/Data/Group fMRI/Group masks/rFG2_r.nii'
%                   };

%% Set files
matlabbatch = {};
matlabbatch{1}.cfg_basicio.cfg_named_file.name = 'Params';
matlabbatch{1}.cfg_basicio.cfg_named_file.files = {{params}};
matlabbatch{2}.cfg_basicio.cfg_named_file.name = 'Space';
matlabbatch{2}.cfg_basicio.cfg_named_file.files = {{spaceImage}};
matlabbatch{3}.cfg_basicio.cfg_named_file.name = 'towarp';
matlabbatch{3}.cfg_basicio.cfg_named_file.files = {images}';

%% Define warp
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1) = cfg_dep;
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).tname = 'Parameter File';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).name = 'class';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).sname = 'Named File Selector: Params(1) - Files';
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.normalise.write.subj.matname(1).src_output = substruct('.','files', '{}',{1});
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep;
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).tname = 'Images to Write';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).name = 'class';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).sname = 'Named File Selector: towarp(1) - Files';
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1).src_output = substruct('.','files', '{}',{1});
matlabbatch{4}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{4}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
                                                          78 76 85];
matlabbatch{4}.spm.spatial.normalise.write.roptions.vox = vx_dim;
matlabbatch{4}.spm.spatial.normalise.write.roptions.interp = 4;
matlabbatch{4}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{4}.spm.spatial.normalise.write.roptions.prefix = 'w';

%% Define reslice
matlabbatch{5}.spm.spatial.coreg.write.ref(1) = cfg_dep;
matlabbatch{5}.spm.spatial.coreg.write.ref(1).tname = 'Image Defining Space';
matlabbatch{5}.spm.spatial.coreg.write.ref(1).tgt_spec{1}(1).name = 'class';
matlabbatch{5}.spm.spatial.coreg.write.ref(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{5}.spm.spatial.coreg.write.ref(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{5}.spm.spatial.coreg.write.ref(1).tgt_spec{1}(2).value = 'e';
matlabbatch{5}.spm.spatial.coreg.write.ref(1).sname = 'Named File Selector: Space(1) - Files';
matlabbatch{5}.spm.spatial.coreg.write.ref(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
matlabbatch{5}.spm.spatial.coreg.write.ref(1).src_output = substruct('.','files', '{}',{1});
matlabbatch{5}.spm.spatial.coreg.write.source(1) = cfg_dep;
matlabbatch{5}.spm.spatial.coreg.write.source(1).tname = 'Images to Reslice';
matlabbatch{5}.spm.spatial.coreg.write.source(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{5}.spm.spatial.coreg.write.source(1).tgt_spec{1}(1).value = 'image';
matlabbatch{5}.spm.spatial.coreg.write.source(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{5}.spm.spatial.coreg.write.source(1).tgt_spec{1}(2).value = 'e';
matlabbatch{5}.spm.spatial.coreg.write.source(1).sname = 'Normalise: Write: Normalised Images (Subj 1)';
matlabbatch{5}.spm.spatial.coreg.write.source(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.spm.spatial.coreg.write.source(1).src_output = substruct('()',{1}, '.','files');

%% Run
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
disp('Warped and resliced.')

end