function [ output_args ] = qc_make_movie_from_nifti(folders,prefix,z,gif)
% Read in all functional volumes corresponding to a session, and then make
% a movie from a  particular slice. Movie will be saved in the quality
% check directory, which will be expected (or created) at the same level as
% the session folders. 
%
% Input: 
% - folders: cell with strings. Each string is the path to a folder that
%            contains the volumes of a session.
% - prefix: string. a wildcard to filter/search volumes with, e.g. '^arf.*\.nii$'
% - z: scalar. slice number. When invalid, will pick a medial one. 
% - gif: binary. If 0, will not export to *.gif. If 1, will save animated
%            gif along movie.
%
% NEEDS: 
% - SPM to read in volumes
%
% WARNING:
% - This function consumes a *lot* of memory, although it clears memory
%   in-between. A hypopriors dataset exceeds 8GB, for example.
%
% Christian Utzerath 2014-15 (Donders Institute)


%% Check or create quality check folder
return_dir = pwd;
cd(folders{1})
cd ..
toplevel = pwd;
qc_dir = [toplevel,filesep,'Quality'];
if ~exist(qc_dir); mkdir(qc_dir); end
cd(return_dir)

%% Read in images, turn into movie
for s = 1:numel(folders)
    session  = folders{s};
    fileseps = strfind(session,filesep);
    sessname = session(fileseps(end)+1:end);
    disp(['Making a slice movie for session: ' sessname ])
    
    clear files volumes Vols dims running_slice slice frame mov 
    
    % Read volumes
    [files] = spm_select('FPListRec',session,prefix);
    volumes = spm_vol(files);
    Vols    = spm_read_vols(volumes);
    
    % Determine target slice
    dims = volumes.dim;
    if z > dims(3); z = ceil(dims(3)/2); end
    
    % Read out target slice from every volume
    running_slice = squeeze(Vols(:,:,z,:));
    
   
    % Convert slices to movie frames. Has to be frame by frame I believe
    clear frame
    gifslice = [];
    for f = 1:size(running_slice,3)
        slice = squeeze(running_slice(:,:,f));
        slice = mat2gray(slice);                % Give grayscale range
        [slice,map] = gray2ind(slice,64);       % Turn into indexed image
        frame(f) = im2frame(slice,map);         % Turn into movie frame        
        gifslice(:,:,1,f) = slice;              % Copy slice into 4D array for animated gif
    end
    
    % Write out gif
    if gif
        imwrite(gifslice,bone,[qc_dir,filesep,sessname,'.gif'],'gif','loopcount',inf,'delaytime',0.05);
        close gcf
    end
    
    % Create and open empty movie
    mov = VideoWriter([qc_dir,filesep,sessname,'.avi']);
    open(mov);
    
    % Write into movie, close file
    writeVideo(mov,frame)    
    close(mov)
    
    disp(['Saved in: ' qc_dir])
end
disp('Done making movies.')

end

