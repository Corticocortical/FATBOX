%% Use Freesurfer to extract V1 masks
% This is an example wrapper that takes you through (most) of the
% steps to create anatomical masks in Freesurfer and export them 
% to a volume (which can be used in SPM). Please note that some
% of the settings and procedures here apply to the environment 
% at the Donders Institute. If you are not a Donderian, things
% might either not work, work differently, or surprisingly work.


%% General information about setting freesurfer up
% You need to make sure that Freesurfer is installed properly on your
% client. This is done through the unix c shell, where you set some global
% environmental variables that Freesurfer will acces. You need to set
% FREESURFER_HOME and SUBJECTS_DIR, although these might also be set
% automatically on your mentat already.
%
% To set, open the unix bash, switch to she c shell (tcsh), and then type:
% setenv FREESURFER_HOME /opt/freesurfer/5.3
% To query that variable, type:
% echo $FREESURFER_HOME
% Note that in c, echo is akin to disp or print. The $ sign states that you
% want to read out the variable that succeeeds the $, so this line reads
% 'Echo [what is stored in] FREESURFER_HOME'. You also need to use the $
% when using variables in other operations. 
%
% Set the subject dir:
% setenv SUBJECTS_DIR /home/predatt/chrutz/ImpexpII/Data/Freesurfer/ANAT
% 
% Verify:
% echo $SUBJECTS_DIR
% echo $FREESURFER_HOME
%
% Then you can tell Freesurfer to set up using those values: 
% source $FREESURFER_HOME/SetUpFreeSurfer.csh
%
% It should return an overview of the most important settings and paths.
%
% Test your installation: 
% freeview -v $SUBJECTS_DIR/S14/mri/orig/001.mgz
%
% Use this script then to convert your subjects' niftis into mgz, copy them
% in the proper folder, and start the procedure with recon-all -subjid S14 -all
%
% Christian Utzerath 2015-16 (Donders Institute)


%% Sample and global settings
sample  = [14:17 19 21 22 24:26 28 30:41 43];
hemis={'l','r'};

copyT1toFS = 0;
labs2vols  = 0;
nativemask = 0;
mnimask    = 1;

selection = 14;
setup_folders
addpath(pipeline_dir)

%% Copy and convert T1 scans into folder structure
if copyT1toFS
disp('Creating Freesurfer-compatible structural images.... ')
for sb = 1:numel(sample)    
    sub = sample(sb);
    selection = sub; setup_folders; 
    disp(['Subject: ' num2str(sub)])
    
    % Get T1
    anatdir = [fMRI_dir,filesep,'Anatomical'];
    cd(anatdir)
    d = dir('s1*.nii')
    T1 = [anatdir,filesep,d.name]
    cd(pipeline_dir)
    
    % Convert
    anat_mgz_file = fullfile(anatdir,'001.mgz');
    cmd_string = sprintf('mri_convert %s %s',T1,anat_mgz_file);
    unix(cmd_string);

    % Copy to Freesurfer/ANAT/Sx
    dest_dir = fullfile(impexp_data_dir,'Freesurfer','subjects',['S' num2str(sub)],'mri','orig');
    if ~exist(dest_dir,'dir')
        mkdir(dest_dir);
    end
    dest = fullfile(dest_dir,'001.mgz');
    cmd_string = sprintf('mv %s %s',anat_mgz_file,dest);
    unix(cmd_string);  
end
disp('===')
end

%% Run reconstruction in c shell
% There should be a script called Reconstruction.csh in the Data/Freesurfer
% folder. Go into the c shell and qsub that script for each subject. 

%% Turn labels into volumes
if labs2vols
labels = {
    'lh.V1.label'
    'rh.V1.label'
    'lh.V2.label'
    'rh.V2.label'
    };

masknames = {
    'V1fs_l.nii'
    'V1fs_r.nii'
    'V2fs_l.nii'
    'V2fs_r.nii'
    };

progress = 0; progressbar(progress);
for s = sample
    % Folder for extracted masks
    labdir = ['/home/predatt/chrutz/ImpexpII/Data/Freesurfer/subjects/S' num2str(s) '/label/extracted'];
    if ~exist(labdir); mkdir(labdir); end;
    
    % Use mri_label2vol to convert labels into masks
    for l = 1:4
        % Which label?
        labelid = ['/home/predatt/chrutz/ImpexpII/Data/Freesurfer/subjects/S' num2str(s) '/label/' labels{l}];
        
        % Regf: space registered to?
        regheader     = ['/home/predatt/chrutz/ImpexpII/Data/Freesurfer/subjects/S' num2str(s) '/mri/aseg.mgz'];
        
        % Seg: segmentation
        segpath  = ['/home/predatt/chrutz/ImpexpII/Data/Freesurfer/subjects/S' num2str(s) '/mri/aseg.mgz'];
        
        % Template for output volume
        anatdir  = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(s) '/fMRI/Anatomical']
        cd(anatdir)
        d = dir('s1*.nii'); T1 = d.name;
        template = [anatdir,filesep,T1];
        
        % Output file
        cd(labdir)
        output = [labdir,filesep,masknames{l}];
        cmd_string = sprintf('mri_label2vol --label %s --temp %s --o %s --regheader %s',labelid,template,output,regheader);
        unix(cmd_string);
    end
    progress = progress +1; progressbar(progress/numel(sample));
end
end

%% MULTIVARAITE MASKS - reslice and make exclusive V1 mask
if nativemask
masknames = {
    'V1fs_l.nii'
    'V1fs_r.nii'
    'V2fs_l.nii'
    'V2fs_r.nii'
    };

% Reslice
progress = 0; progressbar(progress);
for s = sample
    % Define folders
    labdir = ['/home/predatt/chrutz/ImpexpII/Data/Freesurfer/subjects/S' num2str(s) '/label/extracted'];
    anatdir  = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(s) '/fMRI/Anatomical'];
    savedir = [labdir,filesep,'Native']; if ~exist(savedir); mkdir(savedir); end;
    
    % Define images to reslice, and image to reslice to
    spaceImage = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(s) '/fMRI/Results/Final analysis/RSA GLM main task/mask.img'];    
    for m = 1:numel(masknames)
        images_to_reslice{m} = [labdir,filesep,masknames{m}];    
        new_fnames{m} = [labdir,filesep,'r',masknames{m}];
    end
    
    % Reslice
    reslice_images(images_to_reslice,spaceImage)
    
    % Rename and move resliced images to 'Native' folder
    for m = 1:numel(masknames)
        movefile(new_fnames{m},[labdir,filesep,'Native',filesep,masknames{m}])
    end  
   progress = progress +1; progressbar(progress/numel(sample));
end
disp('Native space masks resliced')

% Make a V1 mask that has no V2 voxels
progress = 0; progressbar(progress);
for s = sample
    % Folder
    labdir = ['/home/predatt/chrutz/ImpexpII/Data/Freesurfer/subjects/S' num2str(s) '/label/extracted'];
    savedir = [labdir,filesep,'Native']; 
    targetdir = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(s) '/fMRI/Masks/Native'];
    spaceImage = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(s) '/fMRI/Results/Final analysis/RSA GLM main task/mask.img'];    

    cd(savedir)
    
    % Read V1 and V2 masks
    list={};
    for V = 1:2
        for h = 1:2
            image = [pwd,filesep,['V' num2str(V) 'fs_' hemis{h} '.nii']];
            list{V,h} = make_list_from_image(image);
        end
    end
    
    % Find voxels that are active in both masks
    ind_l = find(list{1,1} == 1 & list{2,1} == 1);
    ind_r = find(list{1,2} == 1 & list{2,2} == 1);
    
    % Remove them from V1 masks
    
    V1fsis_{1} = list{1,1};
    V1fsis_{1}(ind_l) = 0;
    V1fsis_{2} = list{1,2};
    V1fsis_{2}(ind_l) = 0;
    
    % Save V1 isolated mask, copy V1 masks to subject's dir
    for h = 1:2
        fname = [targetdir,filesep,'V1fsis_',hemis{h}];
        dims  = [64 64 26];
        dummy = spaceImage;
        make_binary_mask_from_list(V1fsis_{h},fname,dims,dummy);
    end
    
    % Copy the un-isolated masks    
    for h=1:2
        source = [savedir,filesep,'V1fs_',hemis{h},'.nii'];
        target = [targetdir,filesep,'V1fs_',hemis{h},'.nii'];        
        copyfile(source,target)
    end        
    progress = progress +1; progressbar(progress/numel(sample));
end
end

%% Make MNI space masks to use with multivariate analysis
if mnimask
    progress = 0; progressbar(progress);
    for s = sample
        % Find the directories
        nativedir = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(s) '/fMRI/Masks/Native'];
        tempdir   = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(s) '/fMRI/Masks/MNI/temp'];
        if ~exist(tempdir,'dir'); mkdir(tempdir); end
        MNIdir    = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(s) '/fMRI/Masks/MNI'];
        anatdir   = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(s) '/fMRI/Anatomical'];
        
        
        % Find the masks in native space
        maskpaths = {
            [nativedir,filesep,'V1fs_l.nii']
            [nativedir,filesep,'V1fs_r.nii']
            [nativedir,filesep,'V1fsis_l.nii']
            [nativedir,filesep,'V1fsis_r.nii']
            };
        masknames ={
            'V1fs_l.nii'
            'V1fs_r.nii'
            'V1fsis_l.nii'
            'V1fsis_r.nii'
            };
        
        
        % Copy masks to safe place for processing of the masks
        for i = 1:numel(masknames)
            source = maskpaths{i};
            target = [tempdir,filesep,masknames{i}];
            copyfile(source,target)
        end
        
        % Find parameter file
        clear D
        cd(anatdir)
        D = dir('s*_seg_sn.mat');
        params = [anatdir,filesep,D.name]
        
        % Image in proper space (MNI)
        spaceImage = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(s) '/fMRI/Results/Final analysis/ImpexpII_glm_cat_task/Smoothed/mask.img'];
        
        % Update file paths and get filanames for after processing
        for m = 1:numel(masknames)
            images_to_reslice{m} = [tempdir,filesep,masknames{m}];
            new_fnames{m} = [tempdir,filesep,'rw',masknames{m}];
        end
        warp_and_reslice_images(params,[2 2 2],images_to_reslice,spaceImage)
        
        % Move files upwards
        for m = 1:numel(masknames)
            source = new_fnames{m};
            target = [MNIdir,filesep,masknames{m}];
            copyfile(source,target)
        end
    end
    progress = progress+1; progressbar(progress/numel(sample));
    
end




