function[] = spm_delete_images(folders,remfilter)
% Deletes preprocessed images in a subject's fmri folders.
% 
% Input arguments:
% - folders: will look into these folders (cell string array)
% - remfilter: cell string with one filter per image type to clearn, see
%   below
%
% Useful filters:
% remfilter = {'^f.*\.nii$'};   % remove all converted DICOMs
% remfilter = {'^rf.*\.nii$'};  % remove all realigned files
% remfilter = {'^srf.*\.nii$'}; % remove all smoothed & realigned files
% remfilter = {'.*'};           % remove all files
% C. Utzerath, 2014-15

%% Delete images
setup_folders
fprintf('Removing images for subject %s \n',num2str(selection));
fprintf('============================= \n\n');

matlabbatch = {};
for F=1:length(folders)
    funcdir = folders{F};
    files = {};
    for f=1:numel(remfilter)
        files = [files; cellstr(spm_select('FPList',funcdir,remfilter{f}))];
    end
    
    matlabbatch{F}.cfg_basicio.file_move.files = files;
    matlabbatch{F}.cfg_basicio.file_move.action.delete = false;
end

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);


disp('Deleted all specified files.');