function [beta_vect,beta_files,dims] = read_spm_betas(spmmat_dir,sel_betas,report)
% Turns some beta images into a nvox X ncon matrix.
%
% Input arguments:
% - spmmat_dir: string. Directory with the corresponding SPM.mat (to read info).
% - sel_betas: vector. The numbers of the beta images (eg, 1 for beta_001
%   etc).
% - report: 0/1. Report the names of the conditions (from SPM.mat) or not.
%
% Updated to SPM12 (reads nifti instead of img).
% C. Utzerath, 2014-15

%% Read info from SPM
workdir = pwd;
cd(spmmat_dir)
load SPM

if report==1
    disp(['Reading in the following betas from SPM: ']);
    disp(SPM.xX.name(sel_betas)')
end

dims = SPM.xVol.DIM';

%% Fetch betas
% Find and identify files
%files = spm_select('FPList',pwd,'^beta_.*\.nii$'); % SPM8 command
%files  = cfg_getfile('List',pwd,'^beta_.*\.nii$'); % SPM12 command, gives cell array with files
files  = spm_select('List',pwd,'^beta_.*\.nii$'); % SPM12 command
files  = char(files); % Convert (back) to char array
beta_files = files(sel_betas,:);

% Read out
beta_vol  = spm_vol(files(sel_betas,:));                                % Read in the files corresponding to the selected betas
beta_img  = spm_read_vols(beta_vol);                                    % ""
beta_vect = reshape(beta_img,max(cumprod(dims)),[]);                    % Should be a nvox X 12 matrix
disp('Beta images read.')
cd(workdir)


end

