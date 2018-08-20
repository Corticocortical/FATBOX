function [filenames] = warp_and_reslice(images,deformation,vx_dims,space_image)
% Warp and reslice a set of images between spaces. 
%
% Input argument:
% - images: cell with strings. images to process.
% - deformation: string. Path to (spm12's) deformation field image.
% - vx_dims: 1x3 vector with voxel dimensions of target space
% - space_image: string. path to an image that has the target space. 
%
% Output: 
% - filenames. Cell with path to each image.
%
% Ported for SPM12.

%% Input
% images = {
%     '/home/predatt/chrutz/Hypopriors/Data/P2/MRI/Masks/MNI/temp/rV1_l.nii'
%     '/home/predatt/chrutz/Hypopriors/Data/P2/MRI/Masks/MNI/temp/rV1_r.nii'
%     };
% deformation = '/home/predatt/chrutz/Hypopriors/Data/P2/MRI/Anatomical/iy_s151028081803STD131221107521945416-0010-00001-000192-01.nii';
% vx_dims     = [2 2 2];
% space_image = '/home/predatt/chrutz/Hypopriors/Data/P2/MRI/Results/Paka task GLM/mask.nii';
%
% % C. Utzerath, 2014-17

%% Run
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deformation};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = images;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vx_dims;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{2}.spm.spatial.coreg.write.ref = {space_image};
matlabbatch{2}.spm.spatial.coreg.write.source(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';


spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

%% Return filenames
filenames = {};
for i = 1:numel(images)
    % find slash that marks lowest folder and isolate filename
    sl    = (find(images{i}==filesep)); sl = sl(end);
    fname = images{i}(sl+1:end);
    folder = images{i}(1:sl-1);
    filenames{i} = [folder,filesep,'rw',fname];
end

end

