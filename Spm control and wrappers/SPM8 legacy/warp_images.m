function [ output_args ] = warp_images(parameters,imagestowarp)
% Warps a set of images. Note that depending on what spaces you are working
% with, you might require additional reslicing. Also note that images
% should be properly coregistered beforehand.
%
%
% Input arguments
% - Parameters: cell string {''}. Path to subject's subject-> MNI parameter file
% - images to warp: cell string array {'' ''}. Path to images.
%
% C. Utzerath 2015

%% Deal with input

%% Select files
matlabbatch{1}.cfg_basicio.cfg_named_file.name = 'Parameter';
matlabbatch{1}.cfg_basicio.cfg_named_file.files = {parameters};
matlabbatch{2}.cfg_basicio.cfg_named_file.name = 'ImageToWarp';
matlabbatch{2}.cfg_basicio.cfg_named_file.files = {imagestowarp};

%% Warp                                               
matlabbatch{3}.spm.spatial.normalise.write.subj.matname(1) = cfg_dep;
matlabbatch{3}.spm.spatial.normalise.write.subj.matname(1).tname = 'Parameter File';
matlabbatch{3}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).name = 'class';
matlabbatch{3}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{3}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{3}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).value = 'e';
matlabbatch{3}.spm.spatial.normalise.write.subj.matname(1).sname = 'Named File Selector: Parameter(1) - Files';
matlabbatch{3}.spm.spatial.normalise.write.subj.matname(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.spm.spatial.normalise.write.subj.matname(1).src_output = substruct('.','files', '{}',{1});
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep;
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1).tname = 'Images to Write';
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).name = 'class';
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).value = 'e';
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1).sname = 'Named File Selector: ImageToWarp(1) - Files';
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
matlabbatch{3}.spm.spatial.normalise.write.subj.resample(1).src_output = substruct('.','files', '{}',{1});
matlabbatch{3}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{3}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
                                                          78 76 85];
matlabbatch{3}.spm.spatial.normalise.write.roptions.vox = [2 2 2];
matlabbatch{3}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{3}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.normalise.write.roptions.prefix = 'w';

%% Run
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
disp('Warped.')

end

