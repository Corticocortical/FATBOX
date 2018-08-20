function [ output_args ] = clear_SPM_from_folder(folder)
% This function deletes an SPM.mat file from a specified 
% folder. This is a safety measure in case an SPM.mat is 
% present before an SPM model is run and overwriting is
% not enabled. 
% 
%
% Input:
% - folder: string. folder in which to delete SPM.mat
%
% Christian Utzerath 2014 (Donders Institute)

%% Clear
file = [folder,filesep,'SPM.mat'];    
delete(file) 
disp(['Deleted ' file '.'])


end

