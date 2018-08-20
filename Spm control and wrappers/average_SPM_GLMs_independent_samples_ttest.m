function [] = average_SPM_GLMs_independent_samples_ttest(group,group_names,cons,cn,input,output,mask)
% Use SPM to compare the whole-brain response between two independent
% groups.
%
% Input arguments:
% - group: vecotr. a vector with 0 and 1, describing group membership
% - group_names: cell with strings. names of the two groups.
% - cons: cell string array. File names (not the whole path) of each
%         contrast image that is to be compared.
% - cn:   cell string array. Name for each contrast that is to be compared.
% - input: cell string array. Folders where cons are located.
% - output: string. Where to write results to.
% - mask: explicit mask to be used. string. leave empty if not desired ('')
%
% Compatible with SPM12

%% Batch for group average of this contrast
% % For debugging: give names to input arguments
% group = diagnosis;
% group_names = {'Controls' 'Patients'};
% cons = contrast_files;
% cn = contrast_names;
% input = input_folders;
%
% % C. Utzerath, 2014-16


% Identify members of the two groups 
groups{1} = find(group == 0);
groups{2} = find(group == 1);
    

% Per contrast, assign images per group, output directory, mask, create
% contrasts, and assemble batch file
for icon =  1:numel(cn)        
    % For both groups, find corresponding images for this contrast        
    paths = {};
    paths{1}=cell(numel(groups{1}),1);
    paths{2}=cell(numel(groups{2}),1);
    for g = 1:2        
        for s = 1:numel(groups{g})
            paths{g}{s} = [input{groups{g}(s)},filesep, cons{icon}];
        end
    end
    
    % Assign output directory
    outputdir = [output filesep cn{icon} filesep];
    if ~exist(outputdir); mkdir(outputdir); end
    
    % Test if there is an SPM.mat in the output dir, delete if so
    cd(outputdir)
    spmfile = dir('SPM.mat');
    if isempty(spmfile)
        
    else
        delete(spmfile.name)
    end
    
    % See if mask needed
    if ~isempty(mask)
        mask_entry =[mask,',1'];
        implicitmask = 0;
    else
        mask_entry = '';
        implicitmask = 1;
    end
       
    
    % Create contrast names for both directions of t-test
    contrast_name_1vs2 =  [cn{icon}(1:end-14) group_names{1} ' > ' group_names{2}];
    contrast_name_2vs1 =  [cn{icon}(1:end-14) group_names{2} ' > ' group_names{1}];

    
    matlabbatch = [];
    % Design specification
    matlabbatch{1}.spm.stats.factorial_design.dir = {outputdir};
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = paths{1};
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = paths{2};
    matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;    
    matlabbatch{1}.spm.stats.factorial_design.masking.im = implicitmask;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {mask_entry};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % Estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % Contrasts
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = contrast_name_1vs2;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = contrast_name_2vs1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1;
    

    
    % Run batch for this contrast
    save(['Batch used for ' cn{icon} '.mat'],'matlabbatch')
    spm_jobman('initcfg');         
    spm_jobman('run',matlabbatch);
        
end












% Trashbin
%matlabbatch = [];
% Design specification
%matlabbatch{1}.spm.stats.factorial_design.dir = {outputdir};
%matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = paths;
%     matlabbatch{1}.spm.stats.factorial_design.cov.c = covariate;
%     matlabbatch{1}.spm.stats.factorial_design.cov.cname = covariate_name;
%     matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
%     matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 1;
%     matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
%     matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
%
%     matlabbatch{1}.spm.stats.factorial_design.masking.im = implicitmask;
%     matlabbatch{1}.spm.stats.factorial_design.masking.em = {mask_entry};
%     matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
%     matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
%     matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
%
%     % Estimation
%     matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%     matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
%     matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
%
%     % Contrasts
%     matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%     matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = cn{icon};
%     matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
%     matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.delete = 1;


end