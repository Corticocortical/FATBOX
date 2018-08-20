function [processed_patterns] = rsa_impexpII_pattern_processing(Ylocs,Ytasks,designs,sel_betas,nrun,nregs,fixi_task,fixi_loc,options)
% This function processes response patterns. 
% More or less every step can be toggled on or off, to compare results.
%
% If noise normalization is used (recommended strongly), the rsa toolboxes
% must be on path already for the function to be found.
%
% Input: 
%   - Yloc and Ytasks, scan x voxel matrices from functional scans
%   - designs: paths to SPMs forlocalizer (designs{1}) and main task (designs{2})
%   - sel_betas: vector of indices that tells script which betas to use for
%     from task
%   - nrun: number of sessions in main task
%   - nregs: number of regressors per run
%     for averaging the task data.
%   - fixation_task: vector with indices of all fixation betas
%   - fixation_loc: fixation regressor index for localizer
%   - options: vector of binary toggles for processing options (see below)
%
% Output: 
%   - betapattern: nvox x ncon pattern of beta parameters
%
% Options: [ ... ] 
%       (1): Say what betas you are reading in (determined by SPM)?
%       (2): Use noise normalizatoin (requires rsa toolbox to be ON path)
%       (3): Subtract fixation regressor
%       (4): Z-score the voxels so they have a mean of 0
%       (5): Plot overview of processed data
%    
%
% C. Utzerath, 2014-15

   
%% Deal with input
tell_betas      = options(1);
noise_norm      = options(2);
remove_fix      = options(3);
normalize       = options(4);
plot_processing = options(5);

sel_betamap = reshape(1:numel(sel_betas),nregs,nrun)';

%% Pattern processing
spms{1} = load(designs{1});
spms{2} = load(designs{2});

% Report identities of betas if requested
if tell_betas    
    disp(['Selected betas are: ']);   
    disp(spms{2}.SPM.xX.name(sel_betas)');    
end
Rnames = spms{2}.SPM.xX.name(sel_betas)';

% Use multivariate noise normalization to get best patterns
[u_hat_loc,Sw_hat_loc,resMs_loc,beta_hat_loc] = rsa_noiseNormalizeBeta(Ylocs,spms{1}.SPM);
[u_hat_tas,Sw_hat_tas,resMs_tas,beta_hat_tas] = rsa_noiseNormalizeBeta(Ytasks,spms{2}.SPM);

% Workaround. In one case I have observed that noise normalization produces
% complex values. In those cases, default to beta_hat but produce warning.
if ~isreal(u_hat_loc)
    u_hat_loc = beta_hat_loc;
    %msgbox('Warning: noise normalization of loc betas produced complex values')
    disp('Noise normalization gave complex output, defaulting to betas!')
end
if ~isreal(u_hat_tas)
    u_hat_tas = beta_hat_tas;
    %msgbox('Warning: noise normalization of task betas produced complex values')
    disp('Noise normalization gave complex output, defaulting to betas!')
end



% Pick betas corresponding to conditions of interest
if noise_norm
    loca = u_hat_loc(1:4,:)';
    task = u_hat_tas([sel_betas],:)';
else
    loca = beta_hat_loc(1:4,:)';
    task = beta_hat_tas([sel_betas],:)';
end

% Remove NaNs from betas
patterns = [task loca];
nans = find(isnan(patterns(:,1)));
if numel(nans) > 0
    patterns_nonan = patterns(setdiff(1:size(patterns,1),nans),:);
    disp(['Warning: there are ' num2str(numel(nans)) ' NaNs in your voxel set. Entries removed.']);         
else
    patterns_nonan = patterns;
end

% Subtract fixation if wanted
if remove_fix;
    % Localizer
    fixl = u_hat_loc(fixi_loc,:)';
    loca2 = loca - repmat(fixl,1,size(loca,2));
    %disp(['Subtracted fixation betas from localizer (amplitude: ' num2str(mean(fixl)) ').'])
    
    % Main task
    task2 = task;
    for irun = 1:nrun
        fix = u_hat_tas(fixi_task(irun),:)';
        temp = task(:,sel_betamap(irun,:));
        task2(:,sel_betamap(irun,:)) = temp - repmat(fix,1,size(temp,2));
        %disp(['Subtracted fixation betas from task run (amplitude: ' num2str(mean(fix)) ').'])
    end    
    task = task2;
    loca = loca2;       
else    
    task = task;
    loca = loca;
end
processed_patterns = [task loca];

% Before normalization, average task across runs
taskruns = []; run_inds = reshape([1:size(task,2)],nregs,nrun)';
for irun = 1:4   
    taskruns(:,:,irun) = task(:,run_inds(irun,:));
end
mean_task = mean(taskruns,3); 

% Normalize
if normalize
    patterns_non_norm = [zscore(mean_task,0,2) zscore(loca,0,2)]; % this one produces most 'harmonious' matrices
else
    patterns_non_norm = [mean_task loca];
end
processed_patterns = patterns_non_norm;

% Inspect data processing
if plot_processing
    figure;
    sequentialmap = brewermap(100,'*Blues');
    divergingmap  = brewermap(100,'*Spectral');
    subplot(2,2,1);
    imagesc(patterns_nonan); title('Patterns before normalization'); ylabel('Voxels') ; colormap(sequentialmap);
    subplot(2,2,2);
    imagesc(patterns_non_norm); title('Patterns after normalization & averaging'); ylabel('Voxels'); colormap(sequentialmap);
    subplot(2,2,3);
    imagesc(corr(patterns_nonan,'type','Spearman')); title('Pattern correlation matrix before normalization'); colormap(divergingmap); 
    xlabel('Betas')
    subplot(2,2,4);
    imagesc(corr(patterns_non_norm,'type','Spearman')); title('Pattern correlation matrix after normalization & averaging'); colormap(divergingmap); 
    xlabel('Betas')
    suptitle('Betas and correlations during voxel processing (hot = more)');
end

return

end

