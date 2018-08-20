%% Test out different smoothing parameters with GLMs

%% Test different combinations that create the same FC of ~0.175Hz
%% Set up subject
subjects = {'P1' 'P2' 'P3' 'P4' 'P5' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25' 'S26' 'S27' 'S28' 'S29' 'S30' 'S31' 'S32' 'S33' 'S34' 'S35' 'S36' 'S37' 'S38' 'S39' 'S40' 'S41' 'S42' 'S43' 'S44' 'S45' 'S46' 'S47' 'S48' 'S49' 'S50'};
selected_s = 6;

% Add SPM & installation
installation = '/home/predatt/chrutz/Hypopriors';
addpath(genpath([installation,filesep,'Pipeline']));
spm_dir = add_SPM;

%% Do smoothing and GLM on Paka task
% Settings to compare
combinations = {[11 4] [13 5] [15 6] [19 7]};


for c = 1:numel(combinations)
    % Smooth Sxs's data
    windowSize = combinations{c}(1);
    order      = combinations{c}(2);
    
    folders  = {['/home/predatt/chrutz/Hypopriors/Data/S' num2str(selected_s) '/MRI/EPIs/Run1']
                ['/home/predatt/chrutz/Hypopriors/Data/S' num2str(selected_s) '/MRI/EPIs/Run2']
                ['/home/predatt/chrutz/Hypopriors/Data/S' num2str(selected_s) '/MRI/EPIs/LT localizer']};
    searchstring = ['srf*.nii'];
    dataprop = {[0.68] []}; % TR and task rate
    filprop  = {'savitzky-golay',windowSize,order};
    display = 1;
    fMRI_temporal_smoothing(folders,searchstring,dataprop,filprop,display)
    
    
    % Run GLM on Paka task
    sb = selected_s
    subject = subjects{sb}; t = folder_tree(subject,installation,4);        
    %
    clear glm
    glm.conditions = get_lt_task_regressors(t,0,0);
    glm.batch = [t.pipeline_dir,filesep,'FATBOX/Spm control/SPM batches/generic_GLM_batch.mat'];
    glm.folders = {['/home/predatt/chrutz/Hypopriors/Data/S' num2str(selected_s) '/MRI/EPIs/Run1'] ['/home/predatt/chrutz/Hypopriors/Data/S' num2str(selected_s) '/MRI/EPIs/Run2']};    
    glm.TR = 0.68;
    glm.prefix = '^fsrf.*\.nii$';
    glm.output = ['/home/predatt/chrutz/Hypopriors/Data/S' num2str(selected_s) '/MRI/Results/LT task GLM/Smoothed - SG filter test order ' num2str(order) ' window size ' num2str(windowSize)];
    if ~exist(glm.output); mkdir(glm.output); end
    glm.tcon_names = {'Stim v Fix' 'A v R' 'O v all'};
    glm.tcons      = [-5   1  1  1  1 1;  % Stim v Fix
                       0  -1 -1  1  1 0;  %
                      -1 -1 -1 -1 -1 6];  % ??  
    run_SPM_GLM(glm);
    %
    clear glm
    glm.conditions = get_lt_loca_regressors(t,0,0);
    glm.batch = [t.pipeline_dir,filesep,'FATBOX/Spm control/SPM batches/generic_GLM_batch.mat'];
    glm.folders = {t.lt_localizer_dir};
    glm.TR = 0.68;
    glm.prefix = '^fsrf.*\.nii$';
    glm.output = ['/home/predatt/chrutz/Hypopriors/Data/S' num2str(selected_s) '/MRI/Results/LT loca GLM/Smoothed - SG filter test order ' num2str(order) ' window size ' num2str(windowSize)];
    if ~exist(glm.output); mkdir(glm.output); end
    glm.tcon_names = {'Stim v Fix' 'Obj v Scram'};
    glm.tcons      = [-8 1 1 1 1  1  1  1  1;    % Stim > Fix
                       0 1 1 1 1 -1 -1 -1 -1];   % Obj  > Scrm
    run_SPM_GLM(glm);
    
    
    
    
    
end
