function [] = compare_mask_anat_voxels(group,ivox,contrast,threshold)
%% Put localizer image on anatomy
%
% Displays the masks from the voxel beta analysis and the structural scans,
% each imposed with the subject's localizer t-map. 
%
% Input arguments: 
% - group: 1xn vector of subject numbers, e.g. [5 6 7 8]. N < 15!
% - ivox: voxel selection index from voxel beta analysis. 4 means use 100
%   voxels, 11 means use 500 voxels.
% - contrast: the localizer contrast to draw onto the brains, e.g. 
%   'spmT_0005.img' for objects vs. scrambles.
% - threshold: a value to threshold the localizer image with, e.g. 2. 
%
%
% NEEDS: 
% - Estimated localizer GLM
% - Voxels must already be defined (!! With correct source contrast !!)
% - Access to SPM functions
%
% C. Utzerath, 2014-15


%% Default settings if input not specified
if nargin<4
    if ~exist(group)
        group = [5 6 7 8];
    end
    if ~exist(ivox)
        ivox = 4;  % 4: 100   11: 500
    end
    if ~exist(contrast)
        contrast = 'spmT_0005.img'; % spmT_000t.img -> objects vs. scrambles
    end
    if ~exist(threshold)
        treshold = 2;   % will be divided by 10?
    end
end

%% Break out if too many subjects
if numel(group) > 15
    disp('Registration checks can only be performed on less than 15 subjects at once. Function aborts.')
end

%% Create / find files
% Loop through subjects to create masks. save pathnames for images.
ims = cell(2*numel(group),1);
caps = cell(2*numel(group),1);
tmaps = cell(numel(group),1);

disp('=== Selecting images to display ... ===')
for i = 1:length(group)
   selection = group(i);
   setup_folders
   
   % Create mask, store filename
   create_mask_from_voxel_list(selection,ivox); % create mask from a 100 voxel list   
   ims{i} = [voxel_analysis_folder,filesep,'Localizer mask 100 voxels.nii'];
   caps{i} = ['S' num2str(selection) ' mask'];
   
   % Find anatomical image
   cd([fMRI_dir,filesep,'Anatomical'])
   searchstring = 's*.nii';
   clear search
   search = dir(searchstring);
   
   if isempty(search)
       disp(['Cannot find anatomy for S' num2str(selection) '!']);
       anatfname = uiopen;
   else   
   anatfname = [fMRI_dir,filesep,'Anatomical',filesep,search.name];
   disp(['Found structural: ' anatfname])
   end
   
   cd(pipeline_dir)   
   ims{numel(group) + i}  = anatfname;
   caps{numel(group) + i} = ['S' num2str(selection) ' anatomy']; 
     
   % Find t-maps to display
   tmaps{i} = [cat_localizer_spm_dir,filesep,contrast];
   
end
disp('Masks, structurals, and SPMs selected!')


%% Use check reg and imcalc to display 
images = char(ims);
captions = char(caps);
tmaps = repmat(tmaps,2,1);

% The following code is taken from SPM (slightly modified)
images = spm_vol(images);

spm_figure('GetWin','Graphics');
spm_figure('Clear','Graphics');
spm_orthviews('Reset');

mn = length(images);
n  = round(mn^0.4);
m  = ceil(mn/n);
w  = 1/n;
h  = 1/m;
ds = (w+h)*0.02;

% Plot one image after another
for ij=1:mn
    i = 1-h*(floor((ij-1)/n)+1);
    j = w*rem(ij-1,n);
    handle = spm_orthviews('Image', images(ij),...
        [j+ds/2 i+ds/2 w-ds h-ds]);
    if ij==1, spm_orthviews('Space'); end
    spm_orthviews('AddContext',handle);
    
    % Add a caption
    captions = cellstr(captions);
    mn = numel(captions);
    if ij <= mn
        spm_orthviews('Caption', ij, captions{ij});
    end
    
      
    % Superimpose the localizer results to every image
    % Load contrast image
    sourceimg = tmaps{ij};
    
    % Use imcalc to NaN-out voxels below a threshold
    expression_start = 'i1 + 0./(i1>'
    expression_end   = [num2str(threshold) ')'];
    f = [expression_start expression_end]
    %f = 'i1 + 0./(i1>1)';
    matlabbatch{1}.spm.util.imcalc.input = {sourceimg};
    matlabbatch{1}.spm.util.imcalc.output = 'temp_mask.nii';
    matlabbatch{1}.spm.util.imcalc.outdir = {voxel_analysis_folder};
    matlabbatch{1}.spm.util.imcalc.expression = f;
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    
    % Add that image
    P = [voxel_analysis_folder,filesep,'temp_mask.nii'];
    spm_orthviews('AddColouredImage', ij, P,[1 0 0]);
end
    
end