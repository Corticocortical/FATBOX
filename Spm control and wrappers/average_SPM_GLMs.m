function [] = average_SPM_GLMs(cons,cn,input,output,mask)
% Averages a number of spm (T?) contrasts to look at group responses. 
%
% Input arguments:
% - cons: cell string array. File name (not the whole path) of each contrast
%         image.
% - cn:   cell string array. Name for each contrast.
% - input: cell string array. Folders where cons are located.
% - output: string. Where to write results to.
% - mask: explicit mask to be used. string. leave empty if not desired ('')
%
% Compatible with SPM12
%
% C. Utzerath, 2014-16

%% Batch for group average of this contrast
for icon =  1:numel(cn)    
    % Find files pertaining to this contrast
    paths=cell(numel(input),1);
    for s = 1:numel(input)
        paths{s} = [input{s},filesep, cons{icon}];
    end
   
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
        
    matlabbatch = [];
    % Design specification
    matlabbatch{1}.spm.stats.factorial_design.dir = {outputdir};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = paths;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
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
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = cn{icon};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 1;
    
    % Run batch for this contrast
    spm_jobman('initcfg');         
    spm_jobman('run',matlabbatch);
        
end


end