% This is a pipeline example that wolks through the RSA process. It is
% tailored  to an earlier experiment of ours, so considerable adaptation
% will be required.
%
% This script is able to reproduce the results obtained earlier, when using
% the LOC auto ROI. 
%
% Note: once an ROI has been read in, unless the images or ROI changes, it
% should not be necessary to read it in again. Currently, raw images cannot
% be read from the data file as the data file apparently would exceed
% memory. Think about storing data differently.
%
%
% Resource requirements:
% Takes about 3.5 minutes for one subject and 3 ROIs.
% Note that pattern creation is not optimized for memory. Not sure why
% (because big variables are cleared from memory throughout execution), but
% memory demand can be rather large - use at least 6g if processing whole
% sample.
%
%

%% Initialize
clear all; clc
global analysis_handle;
analysis_handle = 'Final analysis';
addpath(genpath([pwd,filesep,'Helpers']))
tic

%% Set analysis
% What, and in whom, should be looked at?
% Note: In several subjects, RSA noise normalization does not work.
sample  = [14:17 19 21 22 24:26 28 30:41 43];
sample  = 14;
nvoxes = 20:10:500;              % Voxel sets of different size. If an ROI does not have enough voxels, analysis will go through that voxel set only partially.
sample = sample;

% What kind of preparation should be done prior to analysis?
preprocess    = 0;
redo_cons     = 0;
make_patterns = 1;               % Make beta patterns afresh. Needs to be done after changing pattern file. If 0, loads them from the masterfile. Required after changing any ROIs or patterns.
read_scans    = 'yes';            % Read raw data afresh (required for adding or changing ROIs, for example)
analyze       = 1;
report        = 1;

% Set ROIs to process. Suffixes after the ROI name indicate type of ROI 
% and tell the program what file to look for.
% - -no suffix: just read mask from folder
% - -auto: will use the corresponding auto roi contrast to select voxels
%    specify below).
% -  conj: will look for a spherical mask (with -sphere in fname) and limit
%    it by the corresponding contrast specified below. 
ROIs = {'LOC' 'V1fsis'};
hemis = {'l' 'r'};              
prefix = '^arf.*\.nii$';        % Type of image to analyse

% Each ROI will be sorted according to sensitivity in a paricular contrast.
% Per ROI, enter that contrast. 
sens_profile = {'spmT_0003.img' 'spmT_0003.img'}

% Create an evaluable pointer to find mask folder
mask_folder = '[''/home/predatt/chrutz/ImpexpII/Data/S'' num2str(subject) ''/fMRI/Masks/Native'']';

% If you want to run only specific ROIs, set them here
startroi = 1;
endroi   = numel(ROIs);

% Assign measures 
corr_measure = 'spearman';
dist_measure = 'euclidean';

% Design information: where do all the betas go?
nrun = 4; nregs = 16;
active_betas = repmat([ones(1,nregs) zeros(1,20)],1,nrun);  % Betas to read in
sel_betas    = find(active_betas == 1);                     % All selected betas
betamap = reshape(sel_betas,nregs,nrun)';                   % All betas in the experiment
sel_betamap = reshape(1:numel(sel_betas),nregs,nrun)';      % Within (!) selected betas, betas per run
fixi_loc   = 6;              % Index of the fixation regressor in the localizer design
fixi_task  = [17 53 89 125]; % Indices of the fixation regressors in main task

% Where to save data? 
masterfile   = 'RSA_masterfile obj vs scram.mat';
selection = 14; setup_folders;
rsa_store    = [impexp_data_dir,filesep,'Group fMRI',filesep,'RSA'];
rsa_filename = [rsa_store,filesep,masterfile];
figfolder    = [rsa_store,filesep,'Plots'];



%%% Legacy
% For automatic or conjunction ROI masks, enter the source contrasts to
% which the script should refer for each case.
auto_roi_contrasts  = {};
conj_roi_contrasts  = {};

%% Preprocessing, GLM
% if preprocess
%     for subject = sample
%         selection = subject; setup_folders;
%         disp('== Preprocessing localizer for ROI definition.')
%         rsa_run_SPM_preprocessing({},{[fMRI_dir,filesep,'Localizer']},{})
%     end
% end

%% Use this cell to run additional contrasts in the localizer
if redo_cons
    for subject = sample
        selection = subject; setup_folders
        spm_file = {[cat_localizer_spm_dir,filesep,'SPM.mat']};
        
        cnames = {};
        cvecs  = [];
        
        redo_spm_t_contrasts(spm_file,cnames,cvecs);
    end
end

%% Pattern creation (LOC masks look fine, patterns look fine)
if make_patterns
    progressbar('Sample', 'Subject'); done = 0; % pointers for pbar
    for subject = sample
        disp(['Creating patterns for S' num2str(subject) '.'])
        progressbar([],0); % clear ROI/subject counter
        
        % Get this subject prior
        [e u e_tr u_tr pick_mover partners irrelevant] = rsa_impexpII_define_prior(subject);
        
        % Determine folders and SPM designs
        selection = subject; setup_folders;
        session_names = {'Localizer' 'Run1' 'Run2' 'Run3' 'Run4'}; % = folder names
        folders = {}; designs = {};
        for isess = 1:5
            sessdir = [fMRI_dir,filesep,session_names{isess}];
            folders{isess} = sessdir;
        end
        designs{1} = [cat_localizer_spm_dir,filesep,'SPM.mat'];
        designs{2} = [rsa_mainexp_glm_dir,filesep,'SPM.mat'];
        designfolders = {[1]; [2 3 4 5]}; % maps designs to folders
        
        % Read out pattern from mask (task, localizer)
        % Read raw time course; transpose into scans X voxels format foor toolbox
        if strcmp(read_scans,'yes')
            BigYloc = read_volumes_in_roi({folders{designfolders{1}}},prefix,[],'');
            BigYlocs = BigYloc.scans';
            BigYtask = read_volumes_in_roi({folders{designfolders{2}}},prefix,[],'');
            BigYtasks = BigYtask.scans';
            clear BigYloc BigYtask
        end                
        
        % Create patterns, or load them from masterfile
        lists = {}; binary_lists = {};
        for iroi = startroi:endroi
            % Get sensitivity profile for this ROI
            [sens_t] = read_spm_tmap([cat_localizer_spm_dir,filesep,sens_profile{iroi}]);
            [sort_sens,sens_ind] = sort(sens_t,'descend');
            
            disp(['ROI: ' ROIs{iroi} '.'])
            for ihemi = 1:2
                % Determine what mask to use
                myroi = ROIs{iroi};
                if length(myroi) > 5
                    suffix = myroi(end-4   :     end);
                else
                    suffix = '';
                end
                switch suffix
                    case '-auto' % Auto mask from localizer
                        contrastfile = [locafolder,filesep,auto_roi_contrasts{iroi}];
                        [tvec]  = read_spm_tmap(contrastfile);
                        [sorted_tvec,ind_sens] = sort(tvec,'descend');
                        active_voxels_sorted{ihemi} = ind_sens(1:500);
                        
                        % Binarize
                        lists{ihemi} = zeros(length(tvec),1);
                        lists{ihemi}(active_voxels_sorted{ihemi}) = 1;                        
                    case '-conj'
                        % Find mask image
                        maskfile = [conj_roi_masks{iroi} '_' hemis{ihemi} '.nii'];
                        masklist = make_list_from_image(maskfile);
                        
                        % Conjunct with contrast
                        con_list = make_list_from_image([locafolder,filesep,conj_roi_contrasts{iroi}]);                        
                        lists{ihemi} = con_list & masklist;
                    otherwise
                        subjects_mask_dir = eval(mask_folder);
                        lists{ihemi} = make_list_from_image([subjects_mask_dir,filesep,myroi,'_',hemis{ihemi},'.nii']);
                        disp(['Using mask file for ROI: ' myroi ' with ' num2str(length(find(lists{ihemi} > 0))) ' in-mask voxels.']);
                end
                
                % For manual quality control, export result to checkreck im
                export_to_image=0;
                if export_to_image
                    make_binary_mask_from_list(lists{ihemi},[rsa_store,filesep,'roiimage'],[64 64 26],[cat_localizer_spm_dir,filesep,'mask.img'])
                end
                
                % Binarize mask (some masks show p(voxel=grey) by default)
                binary_lists{ihemi} = zeros(length(lists{ihemi}),1);
                ind = find(lists{ihemi} > 0);
                binary_lists{ihemi}(ind) = 1;
                
                % Sort roi by sensitivity profile
                r = ind;                   % thats the roi indices
                s = sens_ind;              % thats the whole brain sorted by sens                                
                w = s(ismember(s,r));      % pos of roi indices in sens profile
                
                % For manual quality check, display sorted
                if export_to_image
                    wl = zeros(numel(sens_ind),1);
                    wl(w) = 1; 
                    make_binary_mask_from_list(wl,[rsa_store,filesep,'roiimage'],[64 64 26],[rsa_loc_glm_dir,filesep,'mask.img'])
                end
                
                % Limit Y to mask, store // Or load from file.
                Yfilename = [rsa_store,filesep,'raw',filesep,'raw_S' num2str(subject) '_' ROIs{iroi} '_' hemis{ihemi} '.mat'];
                if strcmp(read_scans,'yes')
                    Ylocs  = BigYlocs(:,w);
                    Ytasks = BigYtasks(:,w);
                    save(Yfilename,'Ylocs','Ytasks');                                                            
                else                    
                    load(Yfilename)
                end
                
                % Preprocess pattern, then store it
                pattern = rsa_impexpII_pattern_processing(Ylocs,Ytasks,designs,sel_betas,nrun,nregs,fixi_task,fixi_loc,[0 1 1 1 0]);
                clear Ylocs Ytasks
                
                % Load pattern file into working memory, add entry, close
                load(rsa_filename)
                rsa.patterns{subject,iroi,ihemi} = pattern;
                rsa.masks.binary{subject,iroi,ihemi} = binary_lists{ihemi};
                rsa.masks.sorted_active_voxels{subject,iroi,ihemi} = w;
                save(rsa_filename,'rsa');
                disp('Saved beta patterns to rsa file.');
                clear rsa
            end % hemi
            progressbar([],iroi/(1+endroi-startroi)) % update within s counter
        end %roi
        done = done + 1; progressbar(done/numel(sample),[]); % Update progress
        clear  BigYlocs  BigYtasks
    end %subject
end

%% Subject-wise correlation analysis (< 1 minute per ROI)
if analyze
% Load the master file. 
load(rsa_filename);

progressbar('Sample','Subject'); % Init progress bar with 2 counters
done = 0;
for subject = sample
    progressbar([],0); % clear ROI/subject counter
     
    % Get this subject prior for repetitions
    [e u e_tr u_tr pick_mover partners irrelevant] = rsa_impexpII_define_prior(subject);
    
    % For each subject, per stmulus, write down which
    % alternation you want to use for each stimulus template
    % Ie, write down the alternation that 'lands' on the
    % stimulus in question.
    % 1: EA1, 2: EA2, 3: UA1, 4: UA2
    if isEven(subject)
        as = [4;    % corresponds to: bike   -> lion   (ua2)
              2;    % corresponds to: car    -> turtle (ea2)
              3;    % corresponds to: lion   -> bike   (ua1)
              1];   % corresponds to: turtle -> car    (ea1)
        ea  = [2 4];% subject expects stims 2 and 4 to alternate
        ua  = [1 3];
    else
        as = [2;   % corresponds to: bike   -> lion   (ea2)
              4;   % corresponds to: car    -> turtle (ua2)
              1;   % corresponds to: lion   -> bike   (ea1)
              3];  % corresponds to: turtle -> car    (ua1)
        ea  = [1 3];% subject expects stims 1 and 3 to alternate
        ua  = [2 4];
    end
    alt_pick_mover = [0 4 6];
               
    for iroi = startroi:endroi
        for ihemi = 1:2
            % Calculate correlations
            for ivox = 1:numel(nvoxes)
                nvox = nvoxes(ivox);
                if nvox > size(rsa.patterns{subject,iroi,ihemi},1)
                    nvox = size(rsa.patterns{subject,iroi,ihemi},1);
                end
                data      = rsa.patterns{subject,iroi,ihemi}(1:nvox,1:end-4);
                templates = rsa.patterns{subject,iroi,ihemi}(1:nvox,end-3:end);
                
                % Make RDMs / RCMs
                [stim_RDM stim_RCM RCM RDM sc pc sd pd] = deal([]);
                for ipick = 1:3
                    for istim = 1:4
                        % For expected repetitions, choose appropriate pick
                        if ismember(istim,e)
                            rep_beta     = sel_betamap(1,istim+pick_mover(ipick));
                        else
                            rep_beta     = sel_betamap(1,istim);
                        end
                        
                        % Find beta corresponding to relevant alternation
                        % Rotate through partitions if alternation is exp.
                        % I verified correct indices
                        if ismember(istim,ea)
                            this_alt  = 8 + as(istim)+alt_pick_mover(ipick);
                        else
                            this_alt  = 8 + as(istim);
                        end
                        
                        % Define the different patterns
                        R  = data(:,rep_beta);                % Stimulus repetition (task)
                        T  = templates(:,istim);              % Stimulus template (loc.)
                        TU = templates(:,u_tr(istim,2));      % Unexpected stim's template
                        TE = templates(:,e_tr(istim,2));      % Expected stim's template
                        P  = templates(:,partners(istim));    % Partner
                        I1 = templates(:,irrelevant(istim,1));% First irrelevant stim's template
                        I2 = templates(:,irrelevant(istim,2));% Sec. irr. stim's template
                        A  = data(:,this_alt);
                      
                        % Correlate all the different templates
                        RCM(ipick,istim,:,:) = corr([R T P TE TU I1 I2 A],'type',corr_measure);
                        RDM(ipick,istim,:,:) = squareform(pdist([R T P TE TU I1 I2 A]','euclidean'));
                        
                        % Additional measures of interest go here
                        % Does the self/partner-correlation change as function of E?
                        sc(ipick,istim) = RCM(ipick,istim,1,2);
                        pc(ipick,istim) = RCM(ipick,istim,1,3);
                        sd(ipick,istim) = RDM(ipick,istim,1,2);
                        pd(ipick,istim) = RDM(ipick,istim,1,3);
                        
                        % What about distance/correlation to irrelevant stimuli
                        ic(ipick,istim) = mean(RCM(ipick,istim,1,[6 7]));
                        id(ipick,istim) = mean(RDM(ipick,istim,1,[6 7]));
                        
                        % What about the old other-measure
                        oc(ipick,istim) = mean(RCM(ipick,istim,1,[3 6 7]));
                        
                        % Make RCM/RDM that is informative about identities
                        stim_RCM(ipick,istim,:,:) = corr([R templates],'type',corr_measure);
                        stim_RDM(ipick,istim,:,:) = squareform(pdist([R templates]','euclidean'));
                        
                        % Write down correlation between A and T
                        atc(ipick,istim) = RCM(ipick,istim,2,8);
                        atd(ipick,istim) = RDM(ipick,istim,2,8);
                    end % stimulus                                      
                end % pick
                
                % Average stimulus RDMs across picks
                rsa.categories.stim_RDM = squeeze(mean(stim_RDM(:,:,:,:),1));
                rsa.categories.stim_RCM = squeeze(mean(stim_RCM(:,:,:,:),1));
                
                % Describe a category bias (r[ani,ani] > r[ani,vehi]). For
                % each stimulus, show distances to own and other category.
                within_r = []; between_r = [];
                for istim = 1:4
                    if istim < 3 % animals                       
                        within       = setdiff([1 2],istim);
                        between      = [3 4];  
                    elseif istim > 2 % vehicles
                        within       = setdiff([3 4],istim);
                        between      = [1 2];                         
                    end
                    dwithin_r(istim) = rsa.categories.stim_RDM(istim,1,within+1); % +1 because first row/column is correlation of repetition with itself; format is [rep templates(1:4)]
                    dbetween_r(istim) = mean(rsa.categories.stim_RDM(istim,1,between+1));
                    %
                    cwithin_r(istim) = rsa.categories.stim_RCM(istim,1,within+1); % +1 because first row/column is correlation of repetition with itself; format is [rep templates(1:4)]
                    cbetween_r(istim) = mean(rsa.categories.stim_RCM(istim,1,between+1));
                end
                rsa.categories.dwithin(subject,iroi,ihemi,ivox,:)          = dwithin_r;
                rsa.categories.dbetween(subject,iroi,ihemi,ivox,:)         = dbetween_r;
                rsa.categories.dbias(subject,iroi,ihemi,ivox,:)            = dwithin_r - dbetween_r;
                dbias = squeeze(rsa.categories.dbias(subject,iroi,ihemi,ivox,:));
                rsa.categories.dcategory_bias(subject,iroi,ihemi,ivox)   = mean(dbias);
                rsa.categories.dcategory_bias_e(subject,iroi,ihemi,ivox,:) = dbias(e);
                rsa.categories.dcetegory_bias_u(subject,iroi,ihemi,ivox,:) = dbias(u);
                rsa.categories.dbias_modulation(subject,iroi,ihemi,ivox) = mean(dbias(e) - dbias(u));
                %
                rsa.categories.cwithin(subject,iroi,ihemi,ivox,:)          = cwithin_r;
                rsa.categories.cbetween(subject,iroi,ihemi,ivox,:)         = cbetween_r;
                rsa.categories.cbias(subject,iroi,ihemi,ivox,:)            = cwithin_r - cbetween_r;
                cbias = squeeze(rsa.categories.cbias(subject,iroi,ihemi,ivox,:));
                rsa.categories.ccategory_bias(subject,iroi,ihemi,ivox)   = mean(cbias);
                rsa.categories.ccategory_bias_e(subject,iroi,ihemi,ivox,:) = cbias(e);
                rsa.categories.ccetegory_bias_u(subject,iroi,ihemi,ivox,:) = cbias(u);
                rsa.categories.cbias_modulation(subject,iroi,ihemi,ivox) = mean(cbias(e) - cbias(u));
                
                % Average matrices across picks and stimuli.
                RCM_avp   = squeeze(mean(RCM(:,:,:,:),1));   % av. across picks
                RCM_avps  = squeeze(mean(RCM_avp(:,:,:),1)); % av. across stims
                RDM_avp   = squeeze(mean(RDM(:,:,:,:),1));   % av. across picks
                RDM_avps  = squeeze(mean(RDM_avp(:,:,:),1)); % av. across stims
                
                % Store sample-wide vars
                rsa.RCMs.RCM(subject,iroi,ihemi,ivox,:,:) = RCM_avps;
                rsa.RDMs.RDM(subject,iroi,ihemi,ivox,:,:) = RDM_avps;
                
                % Derive more specific correlations of intesrest
                scs = mean(sc,1);
                rsa.RCMs.sce(subject,iroi,ihemi,ivox,:) = mean(scs(e));
                rsa.RCMs.scu(subject,iroi,ihemi,ivox,:) = mean(scs(u));
                %
                pcs = mean(pc,1);
                rsa.RCMs.pce(subject,iroi,ihemi,ivox,:) = mean(pcs(e));
                rsa.RCMs.pcu(subject,iroi,ihemi,ivox,:) = mean(pcs(u));
                %
                ics = mean(ic,1);    
                rsa.RCMs.ice(subject,iroi,ihemi,ivox,:) = mean(ics(e));
                rsa.RCMs.icu(subject,iroi,ihemi,ivox,:) = mean(ics(u));
                
                % How well do alternations correlate to the stimulus templ?
                matc = mean(atc,1);                            % correlation betw. alts and templ.
                mate = mean(atc(ea)); matu = mean(atc(ua));    % correlation betwe. alts and templ., by expectation                                                
                rsa.alt_analysis.matc(subject,iroi,ihemi,ivox,:) = matc;
                rsa.alt_analysis.mate(subject,iroi,ihemi,ivox) = mate;
                rsa.alt_analysis.matu(subject,iroi,ihemi,ivox) = matu;
                rsa.alt_analysis.mat_e_eff(subject,iroi,ihemi,ivox) = mate-matu;
                
                % And how does this compare to self correlations from reps?
                rep_eff  = scs-matc;
                rsa.alt_analysis.rep_effm(subject,iroi,ihemi,ivox) = mean(scs-matc);
                rsa.alt_analysis.rep_eff_exp(subject,iroi,ihemi,ivox) = mean(rep_eff(e) - rep_eff(u));
                
            end % nvox
        end % ihemi
    progressbar([],iroi/(1+endroi-startroi)) % update within s counter
    end % iroi
    
    % Store this subject's data
    rsa.RCMs.info             = 'Dimensions: subject,ROI,hemi,voxels, :';
    disp(['Created RDMs and outcome measures for  S' num2str(subject)'.']);
    done = done + 1; progressbar(done/numel(sample),[]); % Update progress
end % subject

% Save, clean up
save(rsa_filename,'rsa');
toc
end

%% Show, and test, major correlations
if report
load(rsa_filename)

% Which ROI?
pls = []; statsls = [];
for wROI =  startroi:endroi
% Average correlation (across sample) per voxel set
Lscplot = squeeze(nanmean(rsa.RCMs.RCM(sample,wROI,1,:,1,2),1));
Rscplot = squeeze(nanmean(rsa.RCMs.RCM(sample,wROI,2,:,1,2),1));
Lsceplot = squeeze(nanmean(rsa.RCMs.sce(sample,wROI,1,:),1));
Rsceplot = squeeze(nanmean(rsa.RCMs.sce(sample,wROI,2,:),1));
Lscuplot = squeeze(nanmean(rsa.RCMs.scu(sample,wROI,1,:),1));
Rscuplot = squeeze(nanmean(rsa.RCMs.scu(sample,wROI,2,:),1));
%
Licplot = mean(squeeze(nanmean(rsa.RCMs.RCM(sample,wROI,1,:,[6 7]),1)),2);
Ricplot = mean(squeeze(nanmean(rsa.RCMs.RCM(sample,wROI,2,:,[6 7]),1)),2);
Liceplot = squeeze(nanmean(rsa.RCMs.ice(sample,wROI,1,:),1));
Riceplot = squeeze(nanmean(rsa.RCMs.ice(sample,wROI,2,:),1));
Licuplot = squeeze(nanmean(rsa.RCMs.icu(sample,wROI,1,:),1));
Ricuplot = squeeze(nanmean(rsa.RCMs.icu(sample,wROI,2,:),1));
%
Lpcplot = squeeze(nanmean(rsa.RCMs.RCM(sample,wROI,1,:,1,3),1));
Rpcplot = squeeze(nanmean(rsa.RCMs.RCM(sample,wROI,2,:,1,3),1));
Lpceplot = squeeze(nanmean(rsa.RCMs.pce(sample,wROI,1,:),1));
Rpceplot = squeeze(nanmean(rsa.RCMs.pce(sample,wROI,2,:),1));
Lpcuplot = squeeze(nanmean(rsa.RCMs.pcu(sample,wROI,1,:),1));
Rpcuplot = squeeze(nanmean(rsa.RCMs.pcu(sample,wROI,2,:),1));

% Show
overviewplot = figure('Position',[0 0 900 400]);
subplot(1,2,1)
plot(nvoxes,Lscplot,nvoxes,Lsceplot,'b-+',nvoxes,Lscuplot,'b--',nvoxes,Licplot,'r',nvoxes,Liceplot,'r-+',nvoxes,Licuplot,'r--',nvoxes,Lpcplot,'g',nvoxes,Lpceplot,'g-+',nvoxes,Lpcuplot,'g--')
title('Left hemisphere')
legend({'Self' 'Self (e)' 'Self (u)' 'Irr.' 'Irr.(e)' 'Irr(u)' 'Partner' 'Partner(e)' 'Partner(u)'});
xlabel('ROI size (voxels)'); xlabel('Correlation');
cur_ylim = ylim; line([100 100],ylim); ylim(cur_ylim);

subplot(1,2,2)
plot(nvoxes,Rscplot,nvoxes,Rsceplot,'-b+',nvoxes,Rscuplot,'b--',nvoxes,Ricplot,'r',nvoxes,Riceplot,'r-+',nvoxes,Ricuplot,'r--',nvoxes,Rpcplot,'g',nvoxes,Rpceplot,'g-+',nvoxes,Rpcuplot,'g--')
title('Right hemisphere')
legend({'Self' 'Self (e)' 'Self (u)' 'Irr.' 'Irr.(e)' 'Irr(u)' 'Partner' 'Partner(e)' 'Partner(u)'});
xlabel('ROI size (voxels)'); xlabel('Correlation');
tt = ['ROI correlations in ' ROIs{wROI}];
cur_ylim = ylim; line([100 100],ylim); ylim(cur_ylim);
suptitle(tt);


% Save overview plot
sdir = figfolder;
if ~exist(sdir); mkdir(sdir); end
rsa_export_current_figure(sdir,[tt,'.png'])



selection = 14; setup_folders
rmpath(genpath(spm_dir))
for ivox = 1:numel(nvoxes)
    lscetest = rsa.RCMs.sce(sample,wROI,1,ivox);
    rscetest = rsa.RCMs.sce(sample,wROI,2,ivox);
    lscutest = rsa.RCMs.scu(sample,wROI,1,ivox);
    rscutest = rsa.RCMs.scu(sample,wROI,2,ivox);
    [h,pls(wROI,ivox),ci,stattsls(wROI,ivox)] = ttest([lscetest - lscutest],0);
    [h,prs(wROI,ivox),ci,stattsrs(wROI,ivox)] = ttest([rscetest - rscutest],0);
end
addpath(genpath(spm_dir))
end % roi

% Make plot of P's in all ROIs
figure
styles = {'r' 'g' 'b' 'c' 'k' 'r+' 'r-' 'g+' 'g-' 'b+' 'b-' 'ro' 'bo' 'go'};
% left hemi
subplot(1,2,1)
for i = startroi:endroi
    hold on
    plot (nvoxes,pls(i,:),styles{i})
end
hold off
legend(ROIs(startroi:endroi))
title('P-values of expectation effect on self correlations (left hemi)')
xlabel('Amount of voxels included from ROI')
ylabel('P (t-test on difference)')

% right hemit
subplot(1,2,2)
for i = startroi:endroi
    hold on
    plot (nvoxes,prs(i,:),styles{i})
end
hold off
legend(ROIs{startroi:endroi})
title('P-values of expectation effect on self correlations (right hemi)')
xlabel('Amount of voxels included from ROI')
ylabel('P (t-test on difference)')

rsa_export_current_figure(figfolder,['p_per_ROI.png'])
end

