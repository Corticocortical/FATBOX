function [ output_args ] = segment_anatomy(scan)
% Use SPM to segment a subject's stuctural scan.
%
% Input: 
% - scan: string with path to T1.
%
% Takes 2-3 minutes on interactive mentat.
% C. Utzerath, 2014-15

%% Define batch and files
disp('== STARTING SEGMENTATION INTO GREY, WHITE, CSF ==');
matlabbatch = {};
matlabbatch{1}.cfg_basicio.cfg_named_file.name = 'Anatomy';
matlabbatch{1}.cfg_basicio.cfg_named_file.files = {{scan}};

%% Define segmentation
matlabbatch{2}.spm.spatial.preproc.data(1) = cfg_dep;
matlabbatch{2}.spm.spatial.preproc.data(1).tname = 'Data';
matlabbatch{2}.spm.spatial.preproc.data(1).tgt_spec{1}(1).name = 'class';
matlabbatch{2}.spm.spatial.preproc.data(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{2}.spm.spatial.preproc.data(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.spatial.preproc.data(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.spatial.preproc.data(1).sname = 'Named File Selector: Anatomy(1) - Files';
matlabbatch{2}.spm.spatial.preproc.data(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.spatial.preproc.data(1).src_output = substruct('.','files', '{}',{1});
matlabbatch{2}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{2}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{2}.spm.spatial.preproc.output.CSF = [0 0 1];
matlabbatch{2}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{2}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{2}.spm.spatial.preproc.opts.tpm = {
                                               '/home/predatt/chrutz/spm8/tpm/grey.nii'
                                               '/home/predatt/chrutz/spm8/tpm/white.nii'
                                               '/home/predatt/chrutz/spm8/tpm/csf.nii'
                                               };
matlabbatch{2}.spm.spatial.preproc.opts.ngaus = [2
                                                 2
                                                 2
                                                 4];
matlabbatch{2}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{2}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{2}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{2}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{2}.spm.spatial.preproc.opts.biasfwhm = 60;

%% Run
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
disp('Segmented.')

end

