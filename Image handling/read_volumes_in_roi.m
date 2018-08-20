function [roi_scans] = read_volumes_in_roi(folders, prefix,list, filepath)
% When provided with a list of images and a voxel list, this function will
% return an nvox X nscan matrix.
% 
% To read all voxels from the scans, just use a list of all voxels.
%
% Input arguments: 
%       - Folders:  Function will look in each folder for files. Each folder 
%                   is treated as a run, whose scans will be appended to the
%                   output matrix. Cell string array.
%       - prefix:   The prefix of the files to be extracted. Script will
%                   only consider these files. String.
%                   Must have SPM format, e.g. '^srf.*\.nii$'
%       - list:     Nx1 binary vector specifying which voxels are to be
%                   included in the output matrix. If empty, all voxels are
%                   included.
%       - filepath: Full filepath (incl. suffix) to save the output matrix
%                   on hard drive. Leave empty to not save data. Cell
%                   string array.
%
% Output arguments:
%       - roi_scans: a voxel X scans matrix.
%
% C. Utzerath, 2014-15
%
%% Program
disp('Program starts to extract voxel raw data from functional images.')
disp('Reading in images...')
all_funcs = [];

for isess = 1:numel(folders)
    disp(folders{isess});
    
    % Find folders with functional images
    sess_func_dir = folders{isess};
    files = spm_select('FPList',[sess_func_dir],prefix);
    
    % Read in functional images
    func_vol  = spm_vol(files);
    func_img  = spm_read_vols(func_vol);
    
    % Reshape scans into vector format
    dims = func_vol(1).dim;
    nrows = max(cumprod(dims));
    func_vect = reshape(func_img,nrows,size(func_img,4)); % reshape the matrix to a Nvox X Nvol matrix
    funcs{isess} = func_vect;
    disp('Functional images are read in and reshaped.');
    clear func_img;
    
    % Append functional images in one matrix
    all_funcs   = [all_funcs funcs{isess}];
end
disp('Scans collected.')

% If desired, limit all_funcs by list
if ~isempty(list)
    ind = find(list == 1); % slow indexing
    all_funcs = all_funcs(ind,:);
    disp('Scans limited to voxel list.')    
end

% Collect all useful information
roi_scans = struct;
roi_scans.scans  = all_funcs;
roi_scans.ROI    = list;
roi_scans.prefix = prefix;

% If desired, save 
if ~isempty(filepath)
    save(filepath,'roi_scans');
    disp(['Extraced scans saved at: ' filepath]);
end

end
