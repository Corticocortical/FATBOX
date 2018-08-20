function [matlabbatch] = make_batch_7funcs_arf_seg(config)
% Makes an SPM preprocessing batch with the following elements: 
% - Realignment, estimate & reslice
% - Slice timing correction
% - Coregistration T1 to mean functionals
% - Segmentation of T1
%
% This particular batch expects four functional runs. 
%
% Specify the following settings in a config struct variable:
% - func_folders: cell with 1 string per folder path
% - func_prefix:  string with an spm filter to select files (e.g., '^f.*\.nii')
% - T1_dir: string with folder where T1 is stored
% - TR: TR of the EPI sequence
% - n_slices: number of slices
% - slice_order: vector with the slice order for slice time correction
% - reference_slice: reference slice ofr slice time correction
%
%
% Note: with some clever and thorough indexing it should be possible to
%       have this session batches for an arbitrary amount of sessions.

%% Argument mapping
matlabbatch = {};

func_folders = config.func_folders;
prefix = config.func_prefix;
T1_dir = config.T1_dir;
%
n_slices = config.n_slices;
TR = config.TR;
TA = TR - (TR/n_slices);
slice_order = config.slice_order;
reference_slice = config.reference_slice;

%% Files
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {T1_dir};
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^s.*\.nii';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {func_folders{1}};
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^f.*\.nii';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {func_folders{2}};
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^f.*\.nii';
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {func_folders{3}};
matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^f.*\.nii';
matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {func_folders{4}};
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^f.*\.nii';
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {func_folders{5}};
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^f.*\.nii';
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {func_folders{6}};
matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^f.*\.nii';
matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {func_folders{7}};
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^f.*\.nii';
matlabbatch{8}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';

%% Preprocessing
% Realignment and Reslice
matlabbatch{9}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^f.*\.nii)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{9}.spm.spatial.realign.estwrite.data{2}(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^f.*\.nii)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{9}.spm.spatial.realign.estwrite.data{3}(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^f.*\.nii)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{9}.spm.spatial.realign.estwrite.data{4}(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^f.*\.nii)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{9}.spm.spatial.realign.estwrite.data{5}(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^f.*\.nii)', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{9}.spm.spatial.realign.estwrite.data{6}(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^f.*\.nii)', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{9}.spm.spatial.realign.estwrite.data{7}(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^f.*\.nii)', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{9}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{9}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{9}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{9}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{9}.spm.spatial.realign.estwrite.eoptions.interp = 4;
matlabbatch{9}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{9}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{9}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{9}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{9}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{9}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{9}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

% Slice timing
matlabbatch{10}.spm.temporal.st.scans{1}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
matlabbatch{10}.spm.temporal.st.scans{2}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rfiles'));
matlabbatch{10}.spm.temporal.st.scans{3}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 3)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','rfiles'));
matlabbatch{10}.spm.temporal.st.scans{4}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 4)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','rfiles'));
matlabbatch{10}.spm.temporal.st.scans{5}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 5)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','rfiles'));
matlabbatch{10}.spm.temporal.st.scans{6}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 6)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{6}, '.','rfiles'));
matlabbatch{10}.spm.temporal.st.scans{7}(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 7)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{7}, '.','rfiles'));
matlabbatch{10}.spm.temporal.st.nslices = n_slices;
matlabbatch{10}.spm.temporal.st.tr = TR;
matlabbatch{10}.spm.temporal.st.ta = TA;
matlabbatch{10}.spm.temporal.st.so = slice_order;
matlabbatch{10}.spm.temporal.st.refslice = reference_slice;
matlabbatch{10}.spm.temporal.st.prefix = 'a';

% Coreg
matlabbatch{11}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{11}.spm.spatial.coreg.estimate.source(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^s.*\.nii)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{11}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{11}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{11}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{11}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{11}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

% Segment
matlabbatch{12}.spm.spatial.preproc.channel.vols(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^s.*\.nii)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{12}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{12}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{12}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{12}.spm.spatial.preproc.tissue(1).tpm = {'/home/common/matlab/spm12b/tpm/TPM.nii,1'};
matlabbatch{12}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{12}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{12}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{12}.spm.spatial.preproc.tissue(2).tpm = {'/home/common/matlab/spm12b/tpm/TPM.nii,2'};
matlabbatch{12}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{12}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{12}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{12}.spm.spatial.preproc.tissue(3).tpm = {'/home/common/matlab/spm12b/tpm/TPM.nii,3'};
matlabbatch{12}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{12}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{12}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{12}.spm.spatial.preproc.tissue(4).tpm = {'/home/common/matlab/spm12b/tpm/TPM.nii,4'};
matlabbatch{12}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{12}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{12}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{12}.spm.spatial.preproc.tissue(5).tpm = {'/home/common/matlab/spm12b/tpm/TPM.nii,5'};
matlabbatch{12}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{12}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{12}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{12}.spm.spatial.preproc.tissue(6).tpm = {'/home/common/matlab/spm12b/tpm/TPM.nii,6'};
matlabbatch{12}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{12}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{12}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{12}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{12}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{12}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{12}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{12}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{12}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{12}.spm.spatial.preproc.warp.write = [1 1];

end

