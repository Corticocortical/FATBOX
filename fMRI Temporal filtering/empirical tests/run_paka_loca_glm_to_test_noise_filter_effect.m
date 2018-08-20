%% Run S6 Paka loca to test out different noise filters
%Smoothed - moving average lowass 0c2Hz

%%
clear all
subjects = {'P1' 'P2' 'P3' 'P4' 'P5' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25' 'S26' 'S27' 'S28' 'S29' 'S30' 'S31' 'S32' 'S33' 'S34' 'S35' 'S36' 'S37' 'S38' 'S39' 'S40' 'S41' 'S42' 'S43' 'S44' 'S45' 'S46' 'S47' 'S48' 'S49' 'S50'};
selected_s = 6;

% Where dataset is (to be) installed
installation = '/home/predatt/chrutz/Hypopriors';

% Add relevant paths
addpath(genpath([installation,filesep,'Pipeline']))
spm_dir = add_SPM;


%%
sb = selected_s
subject = subjects{sb}; t = folder_tree(subject,installation,4);

% Prepare a GLM (PaKa Localizer, smoothed, for ROI definition?) and run
clear glm
glm.conditions = get_paka_loca_regressors(t,0,0);
glm.batch = [t.pipeline_dir,filesep,'FATBOX/Spm control/SPM batches/generic_GLM_batch.mat'];
glm.folders = {t.paka_localizer_dir};
glm.TR = 0.68;
glm.prefix = '^fsrf.*\.nii$';
glm.output = '/home/predatt/chrutz/Hypopriors/Data/S6/MRI/Results/Paka loca GLM/Smoothed - SG filter approx 0c1hz window 15 order 6';
if ~exist(glm.output); mkdir(glm.output); end
glm.tcon_names = {'Stim v Fix' 'Tri v Ind' 'Ind v Tri'};
glm.tcons      = [-2  1  1;
                   0  1 -1;
                   0 -1  1];
run_SPM_GLM(glm);