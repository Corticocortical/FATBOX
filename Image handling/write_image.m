function [V] = write_image(y, spaceImage, output_filename)
% Write an image to a file, using spm's write vol functions. 
%
%
% Input arguments: 
% - y: either a 3D matrix or an nvox x 1 list with voxel data.
% - space:  full path to an image from with space information are taken.
% - output_filename: full path to where image will be stored.
%
% Example usage: 
% y = ones([64,64,26]);
% spaceImage = [pwd,filesep,'beta_009.nii'];
% output_filename = [pwd,filesep,'written_vol.nii']);
% write_image(v,spaceImage,output_filename);
% 
% C. Utzerath, 2015

%% Deal with input
% Get space info
vol = spm_vol(spaceImage);
vol.dt    = [64 0];
vol.pinfo = [1; 0; 352]; 

% Reshape to matrix if is vector
if numel(size(y)) < 3
    y = reshape(y,vol.dim);
end
    
%% Swap volume information
% Replace filename & description
%k = strfind(output_filename,filesep)
%fname = output_filename(k(end)+1:end);
vol.fname = output_filename;

% Description
vol.descrip = ['Image written out on ' date];

% Write new data in old information
V = spm_write_vol(vol,y);

end

