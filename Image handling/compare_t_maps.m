function [] = compare_t_maps(t_maps,captions)
% Uses SPM's orthview to display a set of t-maps.
% No more than 6 tmaps can be compared at a time for visibility reasons. 
%
% Input arguments: 
% t_maps: cell string array with paths to t-maps
% captions: cell string array with captions for each map. May be left empty
% ({}).
%
% C. Utzerath, 2015

%% Deal with input
if numel(t_maps) > 6
    disp('Program cannot display mor than 6 t maps. Only showing what fits on the screen.')
    t_maps = t_maps{1:6};
end

if isempty(captions) 
    for i = 1:numel(t_maps)
        captions{i} = ['Map ' num2str(i)];
    end
end

%% Program
% Read tmaps
ims = char(t_maps);
images = spm_vol(ims);

% Open canvas, determine formatting
spm_figure('GetWin','Graphics');
spm_figure('Clear','Graphics');
spm_orthviews('Reset');
mn = numel(images);
n  = round(mn^0.4);
m  = ceil(mn/n);
w  = 1/n;
h  = 1/m;
ds = (w+h)*0.02;

% Place items
for ij = 1:mn
    i = 1-h*(floor((ij-1)/n)+1);
    j = w*rem(ij-1,n);
    handle = spm_orthviews('Image', images(ij),...
        [j+ds/2 i+ds/2 w-ds h-ds]);
    if ij==1, spm_orthviews('Space'); end
    spm_orthviews('AddContext',handle);
    
    % Add a caption
    spm_orthviews('Caption', ij, captions{ij});            
end


end
