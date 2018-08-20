function[] = run_SPM_preprocessing(cf)
%% run_spm_preprocessing
% Finds images in the input folders and runs them through spm
% preprocessing.
%
% 
%
% Input arguments:
% - cf: structure with following fields:
%  - subject: name of subject
%  - batch: handle of a function that generates a  matlab batch, e.g. 
%           @make_batch_4funcs_arf_seg
%  - TR: TR of experiment
%  - TA:
%  - slice_order: slice order, for slice time correction
%  - vx_dims: 1x3 vector with voxel dimensions for output scans
%  - FWHM: 1x3 vector with smoothing FWHM
%  - T1_dir: folder where the T1s are to be stored
%  - func_dirs: folders where the functionals should be located, per run*
%     * Raw Dicom files are expeced to reside in subfolders in each folder
%  - flags: 1x2 binary vector with flags to set options: 
%           first flag: perform DICOM conversion
%           second flag: plot subject motion
%           third flag: delete DICOMs after conversion
%
% Example: 
% /fMRI/
% --/T1/ -> T1folder
% ----/DICOM/ -> raw images
% --/Run1/ -> funcfolder1
% ----/DICOM/ -> raw images
% --/Run2/ -> funcfolder2
% ----/DICOM/ -> raw images
%
% cfg.T1_dir    = '/bla/bla/fMRI/T1';
% cfg.func_dirs = {'/bla/bla/fMRI/Run1' '/bla/bla/fMRI/Run2'};
%
% % C. Utzerath, 2014-15

%% General settings and input arguments
skip_scans = 10; % the first 10 scans are wait time and will be discarded (waiting for scanner to reach equilibrium)

% Directories
sessions = {cf.T1_dir cf.func_dirs{:}};
disp(['= Starting Preprocessing  for: ' cf.subject])

% Select preprocessing batch

% Optional flags
if cf.flags(1) == 1
    dicom_import=1;
else
    dicom_import=0;
end

if cf.flags(2) == 1
    plot_head_motion = 1;
else
    plot_head_motion = 0;
end

if cf.flags(3) == 1
    delete_DICOM = 1;
else
    delete_DICOM = 0;
end
    
%% DICOM import
% Files that come out of the scanner have the DICOM format, whereas our
% analysis pipelines use NIFTI. SPM can take care of the import with the
% following batch:

% Create import batch
if dicom_import
disp(['==== DICOM import'])
matlabbatch = {};
for run = 1:numel(sessions)
    raw_dir = [sessions{run},filesep,'DICOM']; % Raw DICOms are stored in child dir of session dir
    outdir  = [sessions{run}];                 % Output in session dir
    
    % Get files for this session
    files = spm_select('FPList',fullfile(raw_dir),'.IMA');                   % list of DICOM files in the input folder
    if run > 1 % Sessions > 1 refer to EPIs, of which the first 4 are skipped
        files(1:skip_scans,:) = [];
    end
    
    % Assemble batch import for this run
    matlabbatch{run}.spm.util.import.dicom.data = cellstr(files);
    matlabbatch{run}.spm.util.import.dicom.root = 'flat';
    matlabbatch{run}.spm.util.import.dicom.outdir = {outdir};
    matlabbatch{run}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{run}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{run}.spm.util.import.dicom.convopts.icedims = 0;    
end
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
disp(['==== DICOMs converted.'])


% Delete DICOMs after conversion
if delete_DICOM
   for run = 1:numel(sessions)
    raw_dir = [sessions{run},filesep,'DICOM']; % Raw DICOms are stored in child dir of session dir
    cd(raw_dir);
    delete([raw_dir,filesep,'*.IMA']);   
   end   
end

end

%% Actual preprocessing (pump files through pre-defined SPM batch)
% This section calls another function to generate a matlab batchjob for
% preprocessing, by handing it the needed settings.
try
    disp(['== Carrying out batch preprocessing for ' cf.subject ' ==']);
    disp(['Batch definition function: ' func2str(cf.batch)])
       
    % Find folders with functional images
    for sess = 2:numel(sessions)        
        funcdir{sess-1} = sessions{sess}; % Functionals to be processed                  
    end
    
    % Folder with structural image(s)    
    anatdir = sessions{1};
    
    % Pass all data needed to the preprocessing batch maker
    cf.preproc.func_folders    = funcdir;
    cf.preproc.T1_dir          = anatdir;
    cf.preproc.func_prefix     = '^f.*\.nii';
    cf.preproc.TR              = cf.TR;
    cf.preproc.n_slices        = cf.n_slices;
    cf.preproc.slice_order     = cf.slice_order;
    cf.preproc.reference_slice = cf.reference_slice;
    cf.preproc.vx_dims         = cf.vx_dims;
    cf.preproc.FWHM            = cf.FWHM;
    
    % Assemble command, run (use a cellfun instead?) to make batch
    command = ['matlabbatch = ',func2str(cf.batch),'(cf.preproc)'];
    eval(command)
    
    % Save a copy of the batch for posterity
%     savedir = [t.func_dir,filesep,'Quality'];
%     if ~exist(savedir); mkdir(savedir); end
%     save([savedir,filesep,cf.subject,'_preprocessing_batch.mat'],'matlabbatch');
%     
    % Run the batch
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    disp(['Preprocessing finished for ', cf.subject '.']);
catch err
    disp('B_preprocessing.m: Something went wrong during preprocessing!');
    rethrow(err)
end

%% Plot subject motion
if plot_head_motion
%     for sess=1:numel(funcdir)
%         func_dir = funcdir{sess};
%         fprintf('Reading: %s ... ',motionfile(1).name);
%         motionfile = dir(fullfile(func_dir,'rp_*.txt'));
%         fid = fopen(fullfile(func_dir,motionfile(1).name),'r');
%         rp = fscanf(fid,' %e %e %e %e %e %e',[6,inf])';
%         fclose(fid);
%         %fprintf('done\n');
%         
%         % plot rp in figure and save
%         hf = figure('Position',get(0,'ScreenSize'));
%         set(hf,'Color','w');
%         
%         subplot(2,1,1);
%         plot(rp(:,1:3));
%         title('Translation','Fontsize',20);
%         ylabel('Distance (mm)','Fontsize',18);
%         xlabel('Time (scans)','Fontsize',18);
%         set(gca,'Fontsize',16);
%         set(gca,'XLim',[1 size(rp,1)]);
%         box off;
%         hl = legend('x-axis','y-axis','z-axis','Location','BestOutside');
%         set(hl,'Fontsize',18);
%         
%         subplot(2,1,2);
%         plot(rp(:,4:6));
%         title('Rotation','Fontsize',20);
%         ylabel('Rotation (deg)','Fontsize',18);
%         xlabel('Time (scans)','Fontsize',18);
%         set(gca,'Fontsize',16);
%         set(gca,'XLim',[1 size(rp,1)]);
%         box off;
%         hl = legend('pitch','roll','yaw','Location','BestOutside');
%         set(hl,'Fontsize',18);
%         
%         figname = fullfile(func_dir,sprintf('realignment_parameters_s%d_%d.jpg',cf.subject,sess));
%         saveas(hf, figname,'jpg');
%     end
end

end % EOF


