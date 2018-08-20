function [ output_args ] = redo_spm_t_contrasts(spm_file, cnames, cvecs)
% Reruns the contrast manager to add or change contrasts. Contrasts are
% replicated across sessions. Contrasts are not tested for validity.
%
% Does only accept t-contrasts at the moment.
%
%
% Input arguments: 
%   - SPM file: full path to the SPM.mat that's to be used
%   - cnames: cell string array with contrast names
%   - cvecs: ncon x nreg matrix for one session (will be replicated)
%
% Generates:
%   - New contrast images in the SPM.mat's directory, potentially
%     overwriting existing ones
%   - New constrasts *should* be added to SPM.mat automatically
%
% Christian Utzerath 2014-15 (Donders Institute)




%% Deal with input arguments
n = numel(cnames);


%% Run contrasts
matlabbatch = {};
% for icon = 1:n
%     matlabbatch{1}.spm.stats.con.spmmat = {spm_file};
%     matlabbatch{1}.spm.stats.con.consess{icon}.tcon.name = cnames{icon};
%     matlabbatch{1}.spm.stats.con.consess{icon}.tcon.convec = cvecs(icon,:);
%     matlabbatch{1}.spm.stats.con.consess{icon}.tcon.sessrep = 'yes';
%     matlabbatch{1}.spm.stats.con.delete = 1;
% end

for icon = 1:n
    matlabbatch{1}.cfg_basicio.cfg_named_file.name = 'SPM';
    matlabbatch{1}.cfg_basicio.cfg_named_file.files = {spm_file};
    
    matlabbatch{2}.spm.stats.con.spmmat(1) = cfg_dep;
    matlabbatch{2}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
    matlabbatch{2}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{2}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{2}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{2}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{2}.spm.stats.con.spmmat(1).sname = 'Named File Selector: SPM(1) - Files';
    matlabbatch{2}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{2}.spm.stats.con.spmmat(1).src_output = substruct('.','files', '{}',{1});
    
    matlabbatch{2}.spm.stats.con.consess{icon}.tcon.name = cnames{icon};
    matlabbatch{2}.spm.stats.con.consess{icon}.tcon.convec = cvecs(icon,:);
    matlabbatch{2}.spm.stats.con.consess{icon}.tcon.sessrep = 'repl';
    matlabbatch{2}.spm.stats.con.delete = 1;
end

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);


end

