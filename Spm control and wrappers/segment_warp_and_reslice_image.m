function [ output_args ] = segment_warp_and_reslice_image(MNI,sT1,image,vx_dims,spaceImage)
%% Warp image into native space, reslice, and copy
% This function warps an image (for example, a mask) into a subject's
% native space, reslices it to have the proper dimensions. The resulting
% file will be placed in the same folder as the to-be-warped-image (copy
% yourself).
%
% Note that the function assumes that the subject's T1 is already in
% register. If you have not run coregistration on the subject's T1 yet, do
% that beforehand.
%
% Function first segments the subject T1 image to obtain normalization and
% inverse normalization parameters. These inverse normalization parameters
% are applied onto the image, with the desired voxel dimensions. The
% resulting image is then resliced to fall into the native subject space,
% using the spaceImage as template. 
%
% Input arguments
% - MNI: cell with 1 string. Path to MNI T1 template file (needed at all?).
% - sT1: cell with 1 string. Path to subject's T1 image.
% - image: cell with 1 string. Path to image to rewarp and reslice.
% - spaceImage: cell with 1 string. Path to image that has the proper
%   space; image will be resliced into this space,
% - vx_dims: 1x3 vector with desired voxel dimensions of warped image.

% C. Utzerath 2015



%% Set files 
disp(['Processing:' image{1}])
matlabbatch = {};
matlabbatch{1}.cfg_basicio.cfg_named_file.name = 'MNI';
matlabbatch{1}.cfg_basicio.cfg_named_file.files = {MNI};
matlabbatch{2}.cfg_basicio.cfg_named_file.name = 'Struct';
matlabbatch{2}.cfg_basicio.cfg_named_file.files = {sT1};
matlabbatch{3}.cfg_basicio.cfg_named_file.name = 'Mask';
matlabbatch{3}.cfg_basicio.cfg_named_file.files = {image};
matlabbatch{4}.cfg_basicio.cfg_named_file.name = 'SpaceDefiningImage';
matlabbatch{4}.cfg_basicio.cfg_named_file.files = {spaceImage};

%% Segment
matlabbatch{5}.spm.spatial.preproc.data(1) = cfg_dep;
matlabbatch{5}.spm.spatial.preproc.data(1).tname = 'Data';
matlabbatch{5}.spm.spatial.preproc.data(1).tgt_spec{1}(1).name = 'class';
matlabbatch{5}.spm.spatial.preproc.data(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{5}.spm.spatial.preproc.data(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{5}.spm.spatial.preproc.data(1).tgt_spec{1}(2).value = 'e';
matlabbatch{5}.spm.spatial.preproc.data(1).sname = 'Named File Selector: Struct(1) - Files';
matlabbatch{5}.spm.spatial.preproc.data(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1});
matlabbatch{5}.spm.spatial.preproc.data(1).src_output = substruct('.','files', '{}',{1});
matlabbatch{5}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{5}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{5}.spm.spatial.preproc.output.CSF = [0 0 0];
matlabbatch{5}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{5}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{5}.spm.spatial.preproc.opts.tpm = {
                                               '/home/common/matlab/spm8/tpm/grey.nii'
                                               '/home/common/matlab/spm8/tpm/white.nii'
                                               '/home/common/matlab/spm8/tpm/csf.nii'
                                               };
matlabbatch{5}.spm.spatial.preproc.opts.ngaus = [2
                                                 2
                                                 2
                                                 4];
matlabbatch{5}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{5}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{5}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{5}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{5}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{5}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{5}.spm.spatial.preproc.opts.msk = {''};

%% Warp
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1) = cfg_dep;
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).tname = 'Parameter File';
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).value = 'e';
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).sname = 'Segment: Norm Params MNI->Subj';
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{6}.spm.spatial.normalise.write.subj.matname(1).src_output = substruct('()',{1}, '.','isnfile', '()',{':'});
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep;
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).tname = 'Images to Write';
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).name = 'class';
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).value = 'e';
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).sname = 'Named File Selector: Mask(1) - Files';
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1).src_output = substruct('.','files', '{}',{1});
matlabbatch{6}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{6}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
                                                          78 76 85];
matlabbatch{6}.spm.spatial.normalise.write.roptions.vox = vx_dims;
matlabbatch{6}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{6}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{6}.spm.spatial.normalise.write.roptions.prefix = 'w';

%% Reslice
matlabbatch{7}.spm.spatial.coreg.write.ref(1) = cfg_dep;
matlabbatch{7}.spm.spatial.coreg.write.ref(1).tname = 'Image Defining Space';
matlabbatch{7}.spm.spatial.coreg.write.ref(1).tgt_spec{1}(1).name = 'class';
matlabbatch{7}.spm.spatial.coreg.write.ref(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{7}.spm.spatial.coreg.write.ref(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{7}.spm.spatial.coreg.write.ref(1).tgt_spec{1}(2).value = 'e';
matlabbatch{7}.spm.spatial.coreg.write.ref(1).sname = 'Named File Selector: SpaceDefiningImage(1) - Files';
matlabbatch{7}.spm.spatial.coreg.write.ref(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1});
matlabbatch{7}.spm.spatial.coreg.write.ref(1).src_output = substruct('.','files', '{}',{1});
matlabbatch{7}.spm.spatial.coreg.write.source(1) = cfg_dep;
matlabbatch{7}.spm.spatial.coreg.write.source(1).tname = 'Images to Reslice';
matlabbatch{7}.spm.spatial.coreg.write.source(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{7}.spm.spatial.coreg.write.source(1).tgt_spec{1}(1).value = 'image';
matlabbatch{7}.spm.spatial.coreg.write.source(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{7}.spm.spatial.coreg.write.source(1).tgt_spec{1}(2).value = 'e';
matlabbatch{7}.spm.spatial.coreg.write.source(1).sname = 'Normalise: Write: Normalised Images (Subj 1)';
matlabbatch{7}.spm.spatial.coreg.write.source(1).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{7}.spm.spatial.coreg.write.source(1).src_output = substruct('()',{1}, '.','files');
matlabbatch{7}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{7}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{7}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{7}.spm.spatial.coreg.write.roptions.prefix = 'r';

%% Run
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

disp(['Done. Look at: rw' image{1}])

end