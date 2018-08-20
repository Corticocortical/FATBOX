%% Searchlight prototype
% Per subject, read in the stimulus beta images, then determine what beta
% reflects what condition, and pass these images into the searchlight in
% order to make a correlation map. 
%
% All RSA GLMs must be performed prior to this (duh).
%
% First searchlight was run with 10 mm sphere, no brain mask.
%
% Note: this file automatically stores results in a folder that is named
%       after the settings used. The group level SPM script has to be set
%       to a folder manually!
%
% C. Utzerath, 2014-15


%% Global ini
clear all; clc
global analysis_handle;
analysis_handle = 'Pilot pipeline 1';
addpath(genpath([pwd,filesep,'Helpers']))
addpath(pwd);

%% Plan analysis
sample  = [14:17 19 21 22 24:26 28 30:41 43];
sampleA = [15:17 19 21 22 24:26 28 30];
sambleB = [31:41 43];
sample = 14;

% What to do?
mode  = 'do';           % do: do the searchlight. read: read results from file (only create and warp images).
use_brain_mask = 1;       % 1: use SPM's brain mask (faster). 0: use all voxels.
radius         = 10;      % Radius of sphere (mm)
warp           = 0;       % Warp output images for further processing
group          = 0;       % Average scans at the end (not working atm due to wonky folder names)

progressbar; progressbar('Subjects...'); progress = 0; progressbar(progress)
for selection = sample;
% Set subject, folder
subject = selection;
setup_folders; addpath(pipeline_dir);

% Create folder based on settings
if use_brain_mask == 1; maskstring = 'with mask'; else maskstring = 'unmasked'; end
foldername = ['Searchlight ' num2str(radius) ' mm ' maskstring];
searchlight_dir = [rsa_dir,filesep,foldername];
if exist(searchlight_dir)
else
    mkdir(searchlight_dir)
end

% Design information: where do all the betas go?
nrun = 4; nregs = 16;
active_betas = repmat([ones(1,nregs) zeros(1,20)],1,nrun);  % Betas to read in
sel_betas    = find(active_betas == 1);                     % All selected betas
betamap = reshape(sel_betas,nregs,nrun)';                   % All betas in the experiment
sel_betamap = reshape(1:numel(sel_betas),nregs,nrun)';      % Within (!) selected betas, betas per run
fixi_loc   = 5;              % Index of the fixation regressor in the localizer design
fixi_task  = [17 53 89 125]; % Indices of the fixation regressors in main task

%% Make beta images
% Define a template volume. Needed for several steps.
example_mask = ([rsa_loc_glm_dir,filesep,'mask.img']);
V = spm_vol(example_mask);          % Keep V as template volume info for later

% Store that mask for later, might use it
maskvol = spm_read_vols(V);
out_of_brain = find(maskvol == 0);

% Define rest of template volume
M = V.mat;                          % Transformation matrix - where voxels go
V.dt = [16 0]; % use float          % Data type set to 32bit floart; SPM beta map default
V.pinfo = [1;0;0];                  % Planar scaling, set to SPM beta map default
V.n     = [1 1];                    % Unknown parameter
V.descrip = 'RSA regressor';        % Fluff


% Create beta images if necessary
if strcmp(mode,'do')    
    cd(rsa_mainexp_glm_dir)
    task_betas = read_spm_betas(pwd,sel_betas);
    cd(rsa_loc_glm_dir)
    loca_betas = read_spm_betas(pwd,[1 2 3 4]);
    
    % Average task betas
    taskruns = [];
    for irun = 1:4
        taskruns(:,:,irun) = task_betas(:,sel_betamap(irun,:));
    end
    mean_task = mean(taskruns,3);
    clear taskruns
    
    % Join datasets (col 1-8: reps; 9-16: alts; 17-20: loca)
    data = [mean_task loca_betas];
    clear task_betas loca_betas taskruns mean_task
    cd(pipeline_dir)
    
    % Write out all the betas for later
    cd(searchlight_dir)
    fs = {};
    for i = 1:size(data,2)
        v  = reshape(data(:,i),V.dim);
        
        % Change the filename and save your own data in that structure
        fs{i} = ['beta_00' num2str(i)  '.nii'];
        V.fname = fs{i};
        spm_write_vol(V,v);
        disp(['Saved image to: ' fs{i}])
    end
end

%% Run searchlight
if strcmp(mode,'do')
    % Assign data for searchlight
    SPM.xY.VY = fs;                                     % Data files (betas)
    if use_brain_mask == 1
        SPM.VM     = example_mask;                          % No mask - all voxels.
    else
        SPM.VM = [];
    end
    searchopt.def  = 'sphere';
    searchopt.spec = [radius];
    fun = 'rsa_impexpII_inspect_searchlight';
    
    % Run searchlight
    [R,debug] = rsa_my_spm_searchlight(SPM,searchopt,fun,subject,sel_betamap);
    save('corr_maps.mat','R')               
    
elseif strcmp(mode,'read')
    cd(searchlight_dir)
    load('corr_maps.mat')
end

% Show self-correlations in preview
figure; for i = 1:26; subplot(5,6,i); imagesc(squeeze(R{1}(:,:,i))); end
suptitle(['Self-correlations (preview S' num2str(selection)])

% Write images (use volume information defined above)
spaceImage = [pwd,filesep,'beta_001.img'];
for ir = 1:numel(R)
    V.fname   = ['corr_' num2str(ir) '.nii']; % File name for image
    V.descrip = 'Correlation map';            % Description
    R{ir}(out_of_brain) = NaN;                % Make out of brain voxels skippable
    VOL = spm_write_vol(V,R{ir});
end

%% Make contrasts
cd(searchlight_dir)

% Define contrasts, and save images
measures = {'SC' 'SCE' 'SCU' 'PC' 'PCE' 'PCU' 'IC' 'ICE' 'ICU'};
cons = {[2 3] [3 2] [5 6] [6 5] [8 9] [9 8]}; % Per con, compute THIS image minus THAT image
c = {}; cfname = {}; 
for icon  = 1:numel(cons)
   c{icon}      =  R{cons{icon}(1)} - R{cons{icon}(2)};
   cfname{icon} = ['con_' measures{cons{icon}(1)} '_vs_' measures{cons{icon}(2)} '.nii'];
   V.fname   = cfname{icon};
   V.descrip = 'Correlation Contrast'; 
   spm_write_vol(V,c{icon});
end
cd(pipeline_dir)

%% Warp contrast and corrs into MNI space
if warp
    % Get segmentation params
    anatdir = [fMRI_dir,filesep,'Anatomical'];
    cd(anatdir)
    d = dir('s*seg_sn.mat'); d = d(1).name;
    pfile = {[anatdir,filesep,d]};
    
    % Get files (cons)
    cd(searchlight_dir)
    images_to_warp = {};
    for  i = 1:numel(cfname)
        images_to_warp{i}  = [searchlight_dir,filesep,cfname{i}];
        
        if i == numel(cfname)
            for j = 1:numel(measures)
                images_to_warp{i+j} = [searchlight_dir,filesep,'corr_',(num2str(j)),'.nii'];
            end
        end
    end
    
    % Warp
    warp_images(pfile,images_to_warp);
    
    
    cd(pipeline_dir)
end

progress = progress+1;
progressbar(progress/numel(sample));
end

%% Try something: average all images  across sample (not working due wo wonky folder names)
if group
sample  = [14 15 16 17 19 21 22 24:26 28 30:41];
path = [impexp_data_dir,filesep,'Group fMRI',filesep,'GroupSearchlight',filesep,'Averaged ' foldername];
if ~exist(path); mkdir(path); end


% Define contrasts and measures again
measures = {'SC' 'SCE' 'SCU' 'PC' 'PCE' 'PCU' 'IC' 'ICE' 'ICU'};
cons = {[2 3] [3 2] [5 6] [6 5] [8 9] [9 8]}; % Per con, compute THIS image minus THAT image

% Collect volumes per subject
[Vcorr Vcon] = deal({});
[VVVCorr VVVCon] = deal([]);
progress = 0; progressbar; progressbar('Collecting files...');  progressbar(progress);
for i = 1:numel(sample)
    selection = sample(i); subject = selection; setup_folders    
    searchlight_dir = [rsa_dir,filesep,'Searchlight'];
    cd(searchlight_dir)        
    
    for icorr = 1:numel(measures)
        Vcorr{icorr} = spm_vol([searchlight_dir,filesep,['wcorr_' num2str(icorr) '.nii']]);
        VVVCorr(:,:,:,i,icorr) = spm_read_vols(Vcorr{icorr});
    end
    
    for icon = 1:numel(cons)
        cfname{icon} = ['con_' measures{cons{icon}(1)} '_vs_' measures{cons{icon}(2)} '.nii'];
        Vs{icon} = spm_vol([searchlight_dir,filesep,['w',cfname{icon}]]);        
        VVVCon(:,:,:,i,icon) = spm_read_vols(Vs{icon});
        cd(pipeline_dir)
    end            
    
    progress = progress + 1; progressbar(progress/numel(sample));    
end
disp('Group data read into WM.')

% Create template volume and folder to store
example_mask = [impexp_data_dir,filesep,'Group fMRI',filesep,'Group Localizer SPM',filesep,'spmT_0001.img']; % Get example space image
V = spm_vol(example_mask);          % Keep V as template volume info for later
V.dt = [16 0];                      % Data type set to 32bit floart; SPM beta map default
V.pinfo = [1;0;0];                  % Planar scaling, set to SPM beta map default
V.n     = [1 1];                    % Unknown parameter
V.descrip = 'Group map';            % Fluff

% Load mask for voxels to ignore
mask = spm_vol([impexp_data_dir,filesep,'Group fMRI',filesep,'Group Localizer SPM',filesep,'mask.img']); % Get example space image
mask = spm_read_vols(mask);
outof_brain = find(mask == 0);

    
% Compute mean images for corrs
cd(path)
for im = 1:numel(measures)
    rmpath(genpath(spm_dir))
    VVVV = squeeze(VVVCorr(:,:,:,:,im));
    VVVm = squeeze(nanmean(VVVV,4));                 % Mean correlation map across subjects
    VVVm(outof_brain) = NaN;                         % Remove voxels outside of brain    
    VVVv = squeeze(nanvar(VVVV,0,4));                 % Variance across subjects
    VVVt = VVVm ./ (VVVv .*  sqrt(numel(sample)));      % Tmap: mean / var*sqrt(N)
    addpath(genpath(spm_dir))
    
    % Writeout
    V.fname = [path,filesep,[measures{im} ' mean.nii']];
    spm_write_vol(V,VVVm);
    V.fname = [path,filesep,[measures{im} ' tmap.nii']];
    spm_write_vol(V,VVVt);
end
    
for ic = 1:numel(cfname)
    rmpath(genpath(spm_dir))
    VVVV = squeeze(VVVCon(:,:,:,:,ic));
    VVVm = squeeze(nanmean(VVVV,4));                 % Mean correlation map across subjects
    VVVm(outof_brain) = NaN;                         % Remove voxels out of brain
    VVVv = squeeze(nanvar(VVVV,0,4));                 % Variance across subjects
    VVVt = VVVm ./ (VVVv .*  sqrt(numel(sample)));      % Tmap: mean / var*sqrt(N)
    addpath(genpath(spm_dir))
    
    % Writeout
    V.fname = [path,filesep,[cfname{ic} ' mean.nii']];
    spm_write_vol(V,VVVm);
    V.fname = [path,filesep,[cfname{ic} ' tmap.nii']];
    spm_write_vol(V,VVVt);
end

    

% Show preview
% figure
% for i = 1: (68-28);
%     slice = i + 17;
%     subplot(5,8,i)
%     imagesc(squeeze(VVVt(:,:,slice)));
%     title(['Slice ' num2str(slice)])
%     colormap('hot')
% end
% suptitle(['Group average of image in ' num2str(numel(sample)) ' subjects'])

end


