% Convert marsbar ROIs to images, lots of them.
%
% % Christian Utzerath 2014-15 (Donders Institute)


sample = [14:17 19 21 22 24:26 28 30:41];

ROI = 'LOC';
hemis = {'l' 'r'};

for subject = sample;
    selection = subject;
    setup_folders
    
    for ihemi = 1:2
        % Location and filename of marsbar ROI object and image for space
        roi_name   = ['S' num2str(selection) '_fROI_' ROI '-blob_' hemis{ihemi}];
        roi_file   = [roi_name '_roi.mat'];
        roi_folder = rsa_dir;
        space_image = [rsa_loc_glm_dir,filesep,'Smoothed',filesep,'mask.img'];
        
        % Load ROI object in workspace
        load([rsa_dir,filesep,roi_file]);
        
        % Save ROI object as image
        marsbar_save_roi_as_image(roi,[rsa_dir,filesep,roi_name,'.nii'],space_image);
        
        disp(['Saved as image: S' num2str(subject)  ', ' roi_name])
    end
end