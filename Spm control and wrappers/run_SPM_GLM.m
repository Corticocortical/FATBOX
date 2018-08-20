function[] = run_SPM_GLM(glm)
%% Run first level model of the category task in SP 
% This is  a wrapper that controls the SPM GLM routine. You can estimate a
% model here, and run t-contrasts on it. F-contrasts are currently not
% supported (compatibility problem from porting from SPM8 to SPM12).
%
% Input arguments: glm, a struct with the following fields:
%  - folders: a cell with strings that are paths to the functional image
%             folders. Motion parameters are expected to be in this folder.
%  - prefix: a spm_select compatible prefix to scan for files within
%             folders, e.g. '^arf.*\.nii$'
%  - TR: TR of the experiment
%  - conditions: per functional run, path to a conditions file
%                that contains the names, onsets, durations variables
%  - output: string, where to save results
%  - tcons: t-contrasts to run (disregard nuisance). ncontrast x
%           ncondition matrix
%  - tcon_names: give names for tcons. 1xncontrast cell with strings
%  - explicit_mask: 'none' or filepath to ane explicit mask to use
%
%% Process input
batch         = glm.batch;
folders       = glm.folders;
TR            = glm.TR;
condfiles     = glm.conditions;
output        = glm.output;
prefix        = glm.prefix;
cvecs         = glm.tcons;
cnames        = glm.tcon_names;
explicit_mask = glm.explicit_mask;

% Process mask
if strcmp(explicit_mask,'none')
    explicit_mask = '';
else explicit_mask = explicit_mask;
end

% Set directory to return to
return_dir = pwd;

%% Assemble GLM batch
load(batch)


% Set output in batch
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(output);

% Walk through all specified sessions/folders, get info, then run
for sess = 1:numel(folders)
    funcdir = folders{sess};
    files = spm_select('FPList',funcdir,prefix);        % Find files that carry the prefix that your preprocessing script saved
    
    % Nuisance regs: load motion, append derivative and squared derivative
    rpfile = spm_select('FPList',funcdir,'^rp.*\.txt$'); % select old movement parameter file
    rp = load(rpfile);
    rp_diff = zeros(size(rp));
    rp_diff(2:end,:) = diff(rp);
    rp_diff2 = rp_diff.^2;
    R = [rp rp_diff rp_diff2];
    filename = 'rp_diff_diff2.mat'; fpath = fullfile(funcdir,'rp_diff_diff2.mat');  save(fpath, 'R');
    rpfile = spm_select('FPList',funcdir,filename); % Load in the new nuisance regressor file
    
    % If an explicit mask is supplied, insert here
    matlabbatch{1}.spm.stats.fmri_spec.mask = {[explicit_mask,',1']};
    
    % Load this session's condition regressors
    condfile = condfiles{sess};
    
    % Set files, nuisance regressors, and condition file in batch
    matlabbatch{1}.spm.stats.fmri_spec.sess(sess).scans = cellstr(files);
    matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi = {condfile};
    matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi_reg = {rpfile};
    
    % Change TR
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    
    % Ripped contrast manager from working batch for SPM12:
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    for c = 1:numel(cnames)
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = cnames{c};
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = cvecs(c,:);
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'repl';        
        matlabbatch{3}.spm.stats.con.delete = 0;
    end         
end


%% Run
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);


%% Proper contrast manager syntax
% matlabbatch{1}.spm.stats.con.spmmat = {'/home/predatt/chrutz/Hypopriors/Data/P1/MRI/Results/Paka task GLM/SPM.mat'};
% matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Stim vs Fix';
% matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [-2 1 1];
% matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
% matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'IC > NoIC';
% matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 1 -1];
% matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'repl';
% matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'IC < NoIC';
% matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 -1 1];
% matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'repl';
% matlabbatch{1}.spm.stats.con.delete = 0;



%% Old F contrast from contrast manager
% Currently, F contrasts do not really get accepted in the contrast manager
% - add them back later. 
%     % Define an overall F contrast
%     load(condfile) % read in condfile to know number of cons
%     overall_F = repmat([diag(ones(1,numel(onsets))) zeros(numel(onsets),size(R,2))],1,4);
%     clear names durations onsets
%     
%     % Enter all contrasts into contrast manager
%     % Contrast 1 is the F contrast. It should always be pre-defined in
%     % the btach, so just overwrite what's written there.
%     matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%     matlabbatch{3}.spm.stats.con.consess{1}.fcon.name =  'Overall F';
%     matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights =  overall_F';
%     matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.delete = 1;
      




end