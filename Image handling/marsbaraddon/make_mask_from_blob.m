%% Make ROI from selected blob
% When pointing at a contrast in the SPM results view, 
% use this short script to write the blob out as an ROI.
%
% Christian Utzerath 2014-15 (Donders Institute)



folder   = rsa_dir;
maskname = ['S' num2str(selection) '_fROI_LOC-blob_l'];
mars_blobs2rois(xSPM, folder, maskname);