%% Easy pattern correlations
%
% This script is not used for RSA per se, but for the related approach of
% simply correlating different voxel patterns. This is less powerful but 
% for exploratory purposes it may be just fine.
%
% This script is a reduced version of my rsa pipeline. I am suspecting that
% all the noise normalization stuff etc. is maybe even hurting me, because
% it's not working 100% reliably anyway. So, here I'm just reading and
% normalizing betas. Rest is identical.
%
% C. Utzerath, 2014-15
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
nvoxes = 10:15:300;              % Voxel sets of different size. If an ROI does not have enough voxels, analysis will go through that voxel set only partially.
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
ROIs = {'LOC' 'V1' 'V2' 'V3'};
hemis = {'l' 'r'};
prefix = '^arf.*\.nii$';        % Type of image to analyse

% Each ROI will be sorted according to sensitivity in a paricular contrast.
% Per ROI, enter that contrast.
sens_profile = {'spmT_0003.img' 'spmT_0003.img' 'spmT_0003.img' 'spmT_0003.img' 'spmT_0003.img' 'spmT_0003.img' }

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
masterfile   = 'RSA_masterfile simple.mat';
selection = 14; setup_folders;
rsa_store    = [impexp_data_dir,filesep,'Group fMRI',filesep,'RSA'];
rsa_filename = [rsa_store,filesep,masterfile];
figfolder    = [rsa_store,filesep,'Plots'];

%%% Legacy
% For automatic or conjunction ROI masks, enter the source contrasts to
% which the script should refer for each case.
auto_roi_contrasts  = {};
conj_roi_contrasts  = {};

%% Read
if make_patterns
    progressbar('Sample', 'Subject'); done = 0; % pointers for pbar
    for subject = sample
        disp(['Creating patterns for S' num2str(subject) '.'])
        progressbar([],0); % clear ROI/subject counter
        
        % Get this subject prior
        [e u e_tr u_tr pick_mover partners irrelevant] = rsa_impexpII_define_prior(subject);
        
        % Read betas
        locafolder = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(subject) '/fMRI/Results/Final analysis/ImpexpII_glm_cat_localizer'];
        taskfolder = ['/home/predatt/chrutz/ImpexpII/Data/S' num2str(subject) '/fMRI/Results/Final analysis/RSA GLM main task'];
        Y_task_betas = read_spm_betas(taskfolder,sel_betas,0);
        Y_task_fixa  = read_spm_betas(taskfolder,[17 53 89 125],0);
        Y_loca_betas = read_spm_betas(locafolder,1:5,0);
        Y_loca_fixa  = read_spm_betas(locafolder,6,0);
        
        % Get the masks
        lists = {}; binary_lists = {};
        for iroi = startroi:endroi
            % Get sensitivity profile for this ROI
            [sens_t] = read_spm_tmap([locafolder,filesep,sens_profile{iroi}]);
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
                
                % Limit betas to mask
                task_betas = Y_task_betas(w,:);
                loca_betas = Y_loca_betas(w,:);
                loca_fixa  = Y_loca_fixa(w,:);
                task_fixa  = Y_task_fixa(w,:);
                
                % Subtract fixation
                % Localizer
                loca2 = loca_betas - repmat(loca_fixa,1,size(loca_betas,2));
                
                % Main task
                task2 = task_betas; task_fixes = [17 53 89 125];
                for irun = 1:4
                    temp = task_betas(:,sel_betamap(irun,:));
                    task2(:,sel_betamap(irun,:)) = temp - repmat(task_fixa(:,irun),1,size(temp,2));
                    %disp(['Subtracted fixation betas from task run (amplitude: ' num2str(mean(fix)) ').'])
                end
                task_betas = task2;
                loca_betas = loca2;
                
                % Before normalization, average task across runs
                taskruns = []; run_inds = sel_betamap;
                for irun = 1:4
                    taskruns(:,:,irun) = task_betas(:,run_inds(irun,:));
                end
                task_betas = mean(taskruns,3);
                
                % Normalize (z-score across columns)
                pattern_norm = [zscore(task_betas,0,2) zscore(loca_betas,0,2)]; % this one produces most 'harmonious' matrices
                pattern_nonorm = [task_betas loca_betas];
                
                % Load pattern file into working memory, add entry, close
                load(rsa_filename)
                rsa.patterns{subject,iroi,ihemi} = pattern_norm;
                rsa.patterns_nonorm{subject,iroi,ihemi} = pattern_nonorm;
                rsa.masks.binary{subject,iroi,ihemi} = binary_lists{ihemi};
                rsa.masks.sorted_active_voxels{subject,iroi,ihemi} = w;
                save(rsa_filename,'rsa');
                disp('Saved beta patterns to rsa file.');
                clear rsa
            end % hemi
            progressbar([],iroi/(1+endroi-startroi)) % update within s counter
        end %roi
        done = done + 1; progressbar(done/numel(sample),[]); % Update progress
    end %subject
end

%% Analyse
% Load the master file.
load(rsa_filename);

progressbar('Sample','Subject'); % Init progress bar with 2 counters
done = 0;
for subject = sample
    progressbar([],0); % clear ROI/subject counter
    % Get this subject prior
    [e u e_tr u_tr pick_mover partners irrelevant] = rsa_impexpII_define_prior(subject);
    for iroi = startroi:endroi
        for ihemi = 1:2
            Y = rsa.patterns{subject,iroi,ihemi};
            inan = find(isnan(Y(:,1)));
            Y(inan,:) = [];
            % Calculate correlations
            for ivox = 1:numel(nvoxes)
                nvox = nvoxes(ivox);
                if nvox > size(Y,1)
                    nvox = size(Y,1);
                end
                
                
                data      = Y(1:nvox,1:end-4);     % reps
                templates = Y(1:nvox,end-4:end-1); % stimulus templates
                scrambles = Y(1:nvox,end);         % scramble pattern
                
                % Make RDMs / RCMs
                [R T TU TE P I1 I2 matrix stim_RDM stim_RCM RCM RDM sc pc sd pd] = deal([]);
                for ipick = 1:3
                    for istim = 1:4
                        % For expected repetitions, choose appropriate pick
                        if ismember(istim,e)
                            rep_beta     = sel_betamap(1,istim+pick_mover(ipick));
                        else
                            rep_beta     = sel_betamap(1,istim);
                        end
                        
                        % Define the different patterns
                        R  = data(:,rep_beta);                % Stimulus repetition (task)
                        T  = templates(:,istim);              % Stimulus template (loc.)
                        TU = templates(:,u_tr(istim,2));      % Unexpected stim's template
                        TE = templates(:,e_tr(istim,2));      % Expected stim's template
                        P  = templates(:,partners(istim));    % Partner
                        I1 = templates(:,irrelevant(istim,1));% First irrelevant stim's template
                        I2 = templates(:,irrelevant(istim,2));% Sec. irr. stim's template
                        matrix = [R T P TE TU I1 I2];
                        
                        % Correlate all the different templates
                        RCM(ipick,istim,:,:) = corr(matrix,'type',corr_measure);
                        RDM(ipick,istim,:,:) = squareform(pdist(matrix','euclidean'));
                        
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
                        
                        % Relate to scrams?
                        dscram(ipick,istim) = pdist([R scrambles]','euclidean');
                        cscr = corr([R scrambles],'type',corr_measure);
                        cscram(ipick,istim) = cscr(1,2);
                        
                        % Make RCM/RDM that is informative about identities
                        stim_RCM(ipick,istim,:,:) = corr([R templates],'type',corr_measure);
                        stim_RDM(ipick,istim,:,:) = squareform(pdist([R templates]','euclidean'));
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
                
                % Derive scramble measures
                c_scr = mean(cscram,1);
                rsa.RCMs.cscre(subject,iroi,ihemi,ivox,:) = mean(c_scr(e));
                rsa.RCMs.cscru(subject,iroi,ihemi,ivox,:) = mean(c_scr(u));
                %
                d_scr = mean(dscram,1);
                rsa.RCMs.dscre(subject,iroi,ihemi,ivox,:) = mean(d_scr(e));
                rsa.RCMs.dscru(subject,iroi,ihemi,ivox,:) = mean(d_scr(u));
                
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

%% Show overview
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
end