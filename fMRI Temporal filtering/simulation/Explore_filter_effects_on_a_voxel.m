%% Explore effects of different temporal filters
% In this script, we will look what different filters exactly do to our
% brain signal on the one hand, and to our design on the other. Look at
% signal pre- and post filtering in both the time and frequency domain.
%
% Goal is to learn the filter properties of the butter and moving average
% filter, and see if we can find a setting that works best a priori.
%
% Usage:
% In the first cell, insert a folder from which functional images can be
% read, and the corresponding searchstring (prefix) for the right images.
% Then specify the TR and the desired properties of the filter(s). The
% script will then read in a voxel time series, a regressor, and show you
% how they look under different filters.
%
% Note that you might want to adjust in the code which voxel and regressor
% are selected.
%
%% Determine data
folder = '/home/predatt/chrutz/Hypopriors/Data/S6/MRI/EPIs/Paka Localizer';
searchstring = ['srf*.nii'];
design = '/home/predatt/chrutz/Hypopriors/Data/S6/MRI/Results/Paka loca GLM/Smoothed - nofilter/SPM.mat';
TR = 0.68;

% Set filter properties
cutoff     = 0.1;       % Cut-off frequency for butterworth low-pass filter
windowSize = 3;         % Window size for moving average filters
k = 2;                  % Polynomial order for the Savitzky-Golay filter
f = 21;                 % Window size for the Savitzky-Golay filter

load(design);
D = SPM.xX.X;





%% Read images, demean, pick voxel
% Identify images (alternatively, use spm_select [ugly file filters])
files = dir([folder,filesep,searchstring]);
filenames = {}; folder_filenames = {};
for f = 1:numel(files);
    filenames{f} = files(f).name;
    folder_filenames{f} = [folder,filesep,files(f).name];
end
slashes = strfind(folder,filesep); session = folder(slashes(end)+1:end); % Write down name of session


% Read in images (note: SPM needs character array for filenames)
disp(['Reading ' num2str(length(files)) ' images in folder: ' folder '...'] )
vol_info  = spm_vol(char(folder_filenames));
volumes   = spm_read_vols(vol_info);
disp('Done.')

% Reshape data into appropriate format for filtering (voxel x time)
vol_dims = size(volumes);
numvox   = max(cumprod(vol_dims(1:3)));
reshaped_volumes = reshape(volumes,numvox,vol_dims(end));

% Demean data (bring every voxel down to zero mean)
volume_means = mean(reshaped_volumes,2);
demeaned_volumes = reshaped_volumes - (repmat(volume_means,1,size(reshaped_volumes,2)));

% Find voxel that is most active, that's your benchmark
contrastimage = read_spm_tmap('/home/predatt/chrutz/Hypopriors/Data/S6/MRI/Results/Paka loca GLM/Smoothed - nofilter/spmT_0001.nii');
[maximum,index] = max(contrastimage); 

% Describe brain data
voxel = demeaned_volumes(index,:);
FS = 1/TR;
T  = 0:TR:(numel(voxel)-1)*TR;

% Describe design data
regressor = D(:,2); % I.e., looking at triangles

% Clear out working memory
clear contrastimage reshaped_volumes demeaned_volumes volumes

%% Get frequency response of unfiltered data and design
% Frequency spectra of unfiltered voxel and regressor
[freqs,spectrum_brain_nofilter] = make_power_density_spectrum(voxel,FS,T,0);
[freqs,spectrum_regressor_nofilter] = make_power_density_spectrum(regressor,FS,T,0);

%% Apply filter: butter
% Filter
[B,A]   = butter(6,0.1); 
voxel_butter  = filter(B,A,voxel);
regressor_butter = filter(B,A,regressor);

% Spectrum
[freqs,spectrum_brain_butter]     = make_power_density_spectrum(voxel_butter,FS,T,0);
[freqs,spectrum_regressor_butter] = make_power_density_spectrum(regressor_butter,FS,T,0);

%% Apply filter: moving average, simple
% Filter voxel
windowSize = 5;                                      % Determine the window size that corresponds best to a particular cutoff frequency
b = (1/windowSize)*ones(1,windowSize); a = 1;        % See formula/example doc filter
voxel_mvavg = filter(b,a,voxel);                     % Filter along columns, row by row (along time)
regressor_mvavg = filter(b,a,regressor);

% Spectrum
[freqs,spectrum_brain_mvavg]     = make_power_density_spectrum(voxel_mvavg,FS,T,0);
[freqs,spectrum_regressor_mvavg] = make_power_density_spectrum(regressor_mvavg,FS,T,0);             

%% Apply filter: moving average, weighted
% Use a Gaussian window
% Filter voxel and regressor with Gaussian window
windowSize = 5; 
b = gausswin(windowSize); a = 1;
voxel_mvavggauss = filter(b,a,voxel);                     % Filter along columns, row by row (along time)
regressor_mvavggauss = filter(b,a,regressor);

[freqs,spectrum_brain_mvavggauss]     = make_power_density_spectrum(voxel_mvavggauss,FS,T,0);
[freqs,spectrum_regressor_mvavggauss] = make_power_density_spectrum(regressor_mvavggauss,FS,T,0);             

%% Apply filter: Savitzky Golai
% This filter operates column by column, so transpose signal

golayvoxel = sgolayfilt(voxel,k,f);
golayreg   = sgolayfilt(regressor,k,f);

[freqs,spectrum_golayvoxel]     = make_power_density_spectrum(golayvoxel,FS,T,0);
[freqs,spectrum_golayreg] = make_power_density_spectrum(golayreg,FS,T,0);             


%% Plotting department
freqlimit = 0.12;


figure
nrows = 5;
ncols = 4;

% Plot unfiltered data, time and frequency
subplot(nrows,ncols,1)
plot(T,voxel); title('Unfiltered voxel signal')
xlabel('Time (s)'); ylabel('Signal amplitude')
subplot(nrows,ncols,2)
plot(T,regressor); title('Unfiltered regressor')
xlabel('Time (s)'); ylabel('Signal amplitude')
%
subplot(nrows,ncols,3)
plot(freqs,spectrum_brain_nofilter); title('Voxel spectrum, unfiltered');
xlim([0 freqlimit]); ylabel('Amplitude'); xlabel('Frequency');
subplot(nrows,ncols,4)
plot(freqs,spectrum_regressor_nofilter); title('Regressor spectrum, unfiltered');
xlim([0 freqlimit]); ylabel('Amplitude'); xlabel('Frequency');

% Plot butter-filtered data
subplot(nrows,ncols,5)
plot(T,voxel_butter); title('Butter-filtered voxel signal')
xlabel('Time (s)'); ylabel('Signal amplitude')
subplot(nrows,ncols,6)
plot(T,regressor_butter); title('Butter-filtered regressor')
xlabel('Time (s)'); ylabel('Signal amplitude')
%
subplot(nrows,ncols,7)
plot(freqs,spectrum_brain_butter); title('Voxel spectrum after butter filter')
xlim([0 freqlimit]); ylabel('Amplitude'); xlabel('Frequency');
subplot(nrows,ncols,8)
plot(freqs,spectrum_regressor_butter); title('Regressor spectrum after butter filter')
xlim([0 freqlimit]); ylabel('Amplitude'); xlabel('Frequency');

% Plot moving-average-filtered data
subplot(nrows,ncols,9)
plot(T,voxel_mvavg)
xlabel('Time (s)'); ylabel('Signal amplitude')
subplot(nrows,ncols,10)
plot(T,regressor_mvavg)
xlabel('Time (s)'); ylabel('Signal amplitude')
%
subplot(nrows,ncols,11)
plot(freqs,spectrum_brain_mvavg); title('Voxel spectrum after moving average filter')
xlim([0 freqlimit]); ylabel('Amplitude'); xlabel('Frequency');
subplot(nrows,ncols,12)
plot(freqs,spectrum_regressor_mvavg); title('Regressor spectrum after moving average filter')
xlim([0 freqlimit]); ylabel('Amplitude'); xlabel('Frequency');

% Plot moving-average-gauss-window-filtered data
subplot(nrows,ncols,13)
plot(T,voxel_mvavggauss)
xlabel('Time (s)'); ylabel('Signal amplitude')
subplot(nrows,ncols,14)
plot(T,regressor_mvavggauss)
xlabel('Time (s)'); ylabel('Signal amplitude')
%
subplot(nrows,ncols,15)
plot(freqs,spectrum_brain_mvavggauss); title('Voxel spectrum after moving average filter (Gaussian window)')
xlim([0 freqlimit]); ylabel('Amplitude'); xlabel('Frequency');
subplot(nrows,ncols,16)
plot(freqs,spectrum_regressor_mvavggauss); title('Regressor spectrum after moving average filter (Gaussian window)')
xlim([0 freqlimit]); ylabel('Amplitude'); xlabel('Frequency');

% Plot golay-filtered data
subplot(nrows,ncols,17)
plot(T,golayvoxel)
xlabel('Time (s)'); ylabel('Signal amplitude')
subplot(nrows,ncols,18)
plot(T,golayreg)
xlabel('Time (s)'); ylabel('Signal amplitude')
%
subplot(nrows,ncols,19)
plot(freqs,spectrum_golayvoxel); title('Voxel after Golay filter')
xlim([0 freqlimit]); ylabel('Amplitude'); xlabel('Frequency');
subplot(nrows,ncols,20)
plot(freqs,spectrum_golayreg); title('Regressor after Golay filter')
xlim([0 freqlimit]); ylabel('Amplitude'); xlabel('Frequency');

suptitle('Frequency spectra')

%% 
