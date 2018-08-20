function [roi_path] = marsbar_make_spherical_roi(centre,radius,roi_label,roi_folder)
% Creates a spherical marsbar ROI object, and saves it.
%
% Input:
% - centre: 3x1 coordinates (mm) for the centre of the ROI
% - radius: scalar value for the radius (mm) of the sphere
% - roi_label: string with the label for the roi. Will also be fileame.
% - roi_folder: directory to store resultant file in.
%
% Output: 
% - roi_path: path to the generated roi, so that roi can be picked up.
%
% Christian Utzerath 2014-15 (Donders Institute)


%% Program
sphere_roi = maroi_sphere(struct('centre',centre,'radius',radius))
saveme     = label(sphere_roi,roi_label);
saveroi(saveme,[roi_folder,filesep,roi_label,'.mat']);


roi_path  = [roi_folder,filesep,roi_label,'.mat'];
end

