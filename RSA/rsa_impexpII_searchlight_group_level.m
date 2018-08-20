% This script is used to run a searchlight.
%
% C. Utzerath, 2014-15


%% Average contrasts using SPM
clear all; clc
sample  = [14:17 19 21 22 24:26 28 30:41 43];
con_names =     {
    'wcorr_1'
    'wcorr_2'
    'wcorr_3'
    'wcorr_4'
    'wcorr_5'
    'wcorr_6'
    'wcorr_7'
    'wcorr_8'
    'wcorr_9'
    'wcon_ICE_vs_ICU'
    'wcon_ICU_vs_ICE'
    'wcon_PCU_vs_PCE'
    'wcon_PCE_vs_PCU'
    'wcon_SCE_vs_SCU'
    'wcon_SCU_vs_SCE'};

%%% ENTER YOUR FOLDERS ON THE SINGLE SUBJECT LEVEL HERE
foldernames{1} = 'Searchlight 10 mm unmasked';
foldernames{2} = 'Searchlight 10 mm with mask';


for f = 2:2
    foldername = foldernames{f};
    group_dir  = ['/home/predatt/chrutz/ImpexpII/Data/Group fMRI/GroupSearchlight/' foldername '/'];
    
    for icon  = 1:numel(con_names)
        gdir = [group_dir,filesep,con_names{icon}];
        mkdir(gdir);
        
        
        % Open batch, find ims
        matlabbatch = {};
        images = {};
        for i = 1:numel(sample)
            sel = sample(i);
            images{i} = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(sel) '/fMRI/Results/Pilot pipeline 1/RSA pipeline 1/RSA/' foldername filesep con_names{icon} '.nii,1'];
        end
        
        % Set images and dirs
        matlabbatch{1}.spm.stats.factorial_design.dir = {[group_dir,filesep,con_names{icon}]};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = images;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        % Estimate
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        
        % Average
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
        matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
        matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = ['Mean ' con_names{icon}];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = 1;
        
        % Run
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);
        
    end % cons
end % folders
