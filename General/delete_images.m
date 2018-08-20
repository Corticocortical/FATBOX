function [ output_args ] = delete_images(folders,filters)
% This function goes through a series of <folders> and deletes a series of
% images defined by <filters> in the filename.
%
% Input: 
% - folders: cell array with 1 string per folder to loop through.
% - filters: cell array with 1 string per prefix to loop through. 
%
% - Example filters:
% 'f*.nii'
% 'arf*.nii'
% 'sarf*.nii'
% 
% Christian Utzerath 2015 (Donders Institute)


%% Proceed
for f = 1:numel(folders)
    cd(folders{f});
    for p = 1:numel(filters)
        pref = filters{p};
        file_list = dir(pref);
        
        delete(pref)
        disp(['Deleted ' num2str(numel(file_list)) ' ' pref ' images in ' folders{f} '.'])
    end
    disp('.')
end


