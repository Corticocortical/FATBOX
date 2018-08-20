function [] = add_t_contrast_to_SPM_GLM(spmmat,cname,cvec)
% This function allows you to add a t-contrast to an existing SPM GLM file.
% This function will call SPM12's batch processing mode to add the
% contrast.
%
% Input argument:
% - spmmat: string. path to the SPM.mat file
% - cname: string. Name of the contrast to add.
% - cvec: string. Vector of the contrast to add. 
%
%
% Requirements:
% - SPM (ideally SPM12) must be on path. 
%
% % C. Utzerath, 2016-17

%% dummy input arguments
spmmat = '/home/predatt/chrutz/Hypopriors/Data/S6/MRI/Results/Paka loca GLM/Smoothed/SPM.mat';
cname  = 'Tri vs Baseline';
cvec   = [-1 1 0 zeros(1,19)];

%% Process
% Load SPM.mat
load(spmmat)

% Get number of existing contrasts
n_cons = numel(SPM.xCon);
cn = n_cons+1;

% Prepare  a batch for a contrast manager 
matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat = {spmmat};
matlabbatch{1}.spm.stats.con.consess{cn}.tcon.name = cname;
matlabbatch{1}.spm.stats.con.consess{cn}.tcon.weights = cvec;
matlabbatch{1}.spm.stats.con.consess{cn}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;

% Run the batch
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);


end

