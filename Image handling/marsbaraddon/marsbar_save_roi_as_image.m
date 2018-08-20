function v = marsbar_save_roi_as_image(o, fname, sp)
% Converts a marsbar ROI object into an ROI image.
%
% Input arguments:
% o    : marsbar ROI object, typically stored and loaded from a mat-file
% fname: complete path and name, including extension, where to store mask
% sp   : image defining the space of the resulting mask image (e.g., tmap=.
% 
% Adapted from Matthew Brett's git hub: https://github.com/matthew-brett

if nargin < 3
      error('Need ROI, output filename and image defining space.');
end

sp = mars_space(sp);
o = maroi_matrix(o, sp);
v = do_write_image(o, fname);