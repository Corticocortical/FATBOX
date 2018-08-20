%% Show normalized voxel values
% This code can be used to show that normalization was succesful. 
% All voxels should have a mean of zero.
%
% C. Utzerath, 2014-15
%

% Open RSA patterns file
fname = [rsa_dir,filesep,'S' num2str(selection) '_' ROI '_patterns.mat']; 
load(fname)

%% Plot conditions
figure; 
for icon = 1:12
subplot(4,3,icon)   
data = rsa.patterns(2).z_betas(:,icon); 
plot(data'); 
title(['Condition ' num2str(icon)])

end

%% Plot voxels
figure
hold on
for i = 1:20
data = rsa.patterns(2).z_betas(i,:); 
plot(data);
title('Condition means for some voxels');
end

%% Plot histogram for a selection of voxels
figure
for i = 1:20
    subplot(5,4,i)
    data = rsa.patterns(2).z_betas(i,:);
    hist(data)
    title(['Histogram of condition values for example voxel ' num2str(i)])
end

%% Plot histogram of each condition
figure
for i = 1:12
    subplot(4,3,i);
    data = rsa.patterns(2).z_betas(:,i);
    hist(data)
    title(['Histogram of normalized voxel values in condition ' num2str(i)])
end


%% Condition and voxel means, histogram of voxel means
rmpath(genpath(spm_dir));
data = rsa.patterns(2).z_betas;

conmeans = nanmean(data,1);
voxmeans = nanmean(data,2);

figure
hist(voxmeans)
title('Histogram of mean values for each voxel in fROI')


addpath(genpath(spm_dir));