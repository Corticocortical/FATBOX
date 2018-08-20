function [] = fMRI_temporal_smoothing(folders,searchstring,dataprop,filprop,display)
%% Temporal filtering/smoothing of fMRI images
% USAGE:
% fMRI_temporal_smoothing(folders,searchstring,dataprop,filprpop,display)
%
% WARNINGS
% 1. Filtering requires a lot of RAM! In my cases, about 900 MB8 images
%    with 2.4mm isometric resolution consume up to 14.5GB RAM during
%    filtering. Reducing this requirement is not easily possible, as you
%    typically have to read in an entire fMRI session at once when relying 
%    on standard toolbox functions (exception could be moving average
%    filters, where it would suffice to read in a few images at a time).
% 2. Mind the dependencies mentioned at the bottom of the documentation.
%
%
% DESCRIPTION
% Filters fMRI datasets temporal to remove high frequencies and/or smooth
% the data temporally. This is particularly useful at high acquisition
% rates (e.g., multiband acceleration), and even more if the task rate is
% slow (such as in blocked designs).
%
% Images produced by this function have an 'f' prefixed to their filename;
% e.g., a NIFTI image that had the prefix 'srf*.nii' would become
% 'fsrf*.nii'.
%
%
% The function offers several different filters to implement:
%   1) Butterworth filter: well-suited for filtering out
%   frequencies, as they have sharp roll-offs and good stopband attenuation
%   in the frequency domain. They are less suited for smoothing of fMRI
%   signals, as they tend to shift BOLD peaks and distort their shape.
%   For these low-pass filters, an order and cutoff frequency must be
%   supplied.
%
%   2) Savitzky-Golay filter: well-suited for smoothing the BOLD response
%   temporally, as the filter's symmetrical impulse repsonse preserves the
%   BOLD peak shapes. The frequency response is less optimal: while it has
%   near-perfect pass-bands and good roll-off, stopband attenuation is
%   mediocre. The filter differs from moving average filters in that a
%   polynomial is fitted to each window segment of a time series (if the
%   order of the polynomial is 1, the filter is identical to a moving
%   average filter). The choice of order and window size comes with
%   trade-offs: larger window sizes result in smoother signal, but
%   attenuate the entire time series (i.e. signal loss). Higher-order
%   polynomials can fit more complex waveforms, but are consequently less
%   effective at removing high frequency noise.
%   Window size and polynomial order must be supplied.
%
%   Read Schaefer (2011). What is a Savitzky-Golay filter? IEEE Signal
%   processing magazine, 28(4) for a brief overview of this filter.
%
%   3) Moving average filter: simply slides a window along the time series,
%   within which all data points are averaged. While this can lead to
%   smooth signals, the 'boxy' nature of this filter creates hefty ripples
%   in the frequency domain of the time series. Window size must be
%   supplied.
%
%   4) Moving average Gaussian window filter: functions like a moving
%   average filter, but instead of a box car, a Gaussian is sled over the
%   time series, i.e. every data point gets a weight assigned from a
%   Gaussian distribution that is centered about the middle of the current
%   window. Improved frequency response as opposed to standard moving
%   average filter.
%
%
% Input arguments:
% - folders: cell array containing strings. Each string is the path of a
%            folder with images. Function will take all images therein.
% 
% - search_string: string. Function will read in images that satisfy the 
%            search string. Format same as in MATLAB/Linux bash. 
%            Example: 'srf*.nii'.
% 
% - dataprop: cell array containing, in each cell, properties of the data.
%       dataprop{1}: scalar. TR of the experiment in seconds (e.g., 0.68).
%       dataprop{2}: (optional) scalar. Task rate. Only meaningful if an
%                    unjittered design is used. If supplied, will be used
%                    to mark events on the output plot.
%
% - filprop: cell array containing, in each cell, properties of desired
%            filter.
%       filprop{1}: string. Filter type. Can be 'Butterworth', 'Moving
%                   average', 'Gaussian Moving Average', 'Savitzky-Golay'.
%                   Not case sensitive. 
%       filprop{2}: First parameter. Depending on filter type:
%                   If Butterworth: cut-off frequency
%                   If moving average: window size
%                   If gaussian moving average: window size
%                   If Savitzky-Golay: window size
%       filprop{3}: Second parameter. Depending on filter type:
%                   If Butterworth: filter order
%                   If moving average: unused, skip or leave empty
%                   If Gaussian mov. aver.: weights per sample within window
%                   If Savitzky-Golay: order of fitted polynomial
%
% - display: binary. If 1, will output and save a comparison of filtered
%            and unfiltered data, if 0, not. If task rate is supplied, will
%            mark events on plots.
%
% EXAMPLES
% fMRI_temporal_smoothing({'run1/images/' 'run2/images/'},'srf*.nii',{[0.68] [14]},{'Savitzky-Golay',11,3},1)
% 
% This will smooth all NIFTI images whose names start with 'srf' and are
% found in the folders of run 1 and 2. Images are supposed to have a TR of
% 0.68 and will be smoothed using a Savitzky-Golay filter with a polynomial
% order of 3 and a window size of 11. A comparison of the filtered and
% unfiltered timeseries with an example voxel will be displayed and saved
% along with the smoothed images. Smoothed images are saved in the same
% folder as the original images.
%
%
% DEPENDENCIES
% - SPM: to read and write images, this function uses SPM. Tested with
%        SPM12, but should be compatible with other versions as well.
% - DSP toolbox: The filtering is done through functions that are in 
%        Matlab's signal processing toolbox. If unavailable, exchange
%        through equivalent functions.
% - progressbar: Function refers twice to an external function called
%                progressbar to display a prettier wait bar than matlab's
%                default one. If you do not want to use this function, you
%                can safely comment/delete these calls.
%
%
% RELATED FUNCTIONS
% Several other functions in this folder can be, or are used, to determine
% for example what window size might fit to a desired cutoff frequency, to
% plot spectra, and so on. All functions should be documented.
%
% chrutz 2016
%                               
%

%% Create input arguments
% clear all; clc
%
% Generalized input argument(s)
% folders  = {'/home/predatt/chrutz/Hypopriors/Data/S6/MRI/EPIs/Paka Localizer'};
% searchstring = ['srf*.nii'];
% dataprop = {[0.68] [13.6014]}; % TR and task rate
% filprop  = {'savitzky-golay',15,6};
% display = 1;

%% Check if SPM/Fieldtrip are on path and whether filtering makes sense
disp('== FMRI TIME SERIES TEMPORAL FILTERING ==')
tic_call = tic;
disp('==========================================================')
if ~exist('spm')
    disp('Could not find SPM on path. Trying to add SPM using a function...')
    try
        add_SPM;
        disp('Used add_SPM function to add SPM to path')
    catch
        disp('Could not add SPM using add_SPM function - does it exist?')
    end    
else
    disp('Checked for SPM and found it.')
end
disp('   Errors when reading volumes can stem from incorrect initialization of SPM.')

% Signal properties
TR = dataprop{1};        % TR
if ~isempty(dataprop{2})
    task_rate = dataprop{2}; % Task rate (optional)
    task_frequency = 1/task_rate;              
end
FS      = 1/TR;          % Sampling frequency in Hertz
Nyquist = FS/2;          % Nyquist frequency
  

% Filter properties, depending on filter choice
switch lower(filprop{1})
    case lower('butterworth')
        disp('Filtering with butterworth filter.')
        ftype  = 'butter';
        cutoff = filprop{2}/FS; % the filer expects the cutoff to be supplied in 0<Fc<1, whereby 1 == FS
        order  = filprop{3};
    case lower('moving average')
        disp('Filtering with moving average.')
        disp('   Use moving_average_window_size_lowpass.m to find window sizes that approximately fit certain frequency cutoffs.');
        ftype = 'mvavg';
        windowSize = filprop{2};
    case lower('gaussian moving average')
        disp('Filtering with moving average in Gaussian window.')
        disp('   Use moving_average_window_size_lowpass.m to find window sizes that approximately fit certain frequency cutoffs.');
        ftype = 'mvavggauss';
        windowSize = filprop{2};
        weights    = filprop{3};
    case lower('savitzky-golay')        
        disp('Filtering with Savitzky-Golay filter.')
        ftype = 'golay';
        windowSize = filprop{2};
        order      = filprop{3};
        if isEven(windowSize)
            windowSize = windowSize + 1;
            disp('   Window size was increased by 1 to make it uneven (see doc for Savitzky-Golay filters).')
        end
    otherwise
        disp('Unknown filter selected! Cannot proceed.')
        return
end


% To estimate a periodic signal you need twice as many samples as Hertz
% (up, down, up, down)
if strcmp(ftype,'butter')
    if filprop{2} > Nyquist
        disp('Warning: the cutoff frequency is higher than the dataset''s Nyquist frequency. Filtering may fail!')
    elseif filprop{2} < Nyquist
        disp(['   Nyquist frequency of data is ' num2str(Nyquist) 'Hertz; chosen cutoff is below that at ' num2str(filprop{2}) ' Hertz.'])
    else
    
    end
end


%% Read in a series of images, filter, save 
for F = 1:numel(folders)
    disp('.')
    filtertic = tic;
    % Identify images (alternatively, use spm_select [ugly file filters])
    folder = folders{F};
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
    volumes = reshape(volumes,numvox,vol_dims(end));
    reshaped_vol_dims = size(volumes);

    % If display is desired, save one voxel for display purposes
    if display; voxel = 1001; pre_filter_signal = volumes(voxel,:); end
    
    % Demean data (bring every voxel down to zero mean)
    volume_means = mean(volumes,2);
    volumes = volumes - (repmat(volume_means,1,reshaped_vol_dims(2)));
   
    % Apply filter
    if strcmp(ftype,'mvavg')
        % Perform temporal smoothing by computing a moving average of the
        % signal, using the matlab example for filter function (appr. 28s)
        disp('Performing temporal smoothing via moving average.')
        b = (1/windowSize)*ones(1,windowSize); a = 1;       % See formula/example doc filter
        volumes = filter(b,a,volumes,[],2);        % Filter along columns, row by row (along time)
    elseif strcmp(ftype,'butter')
        % Otherwise just filter the data as one time series (approx. 10s)
        disp(['Running Butterworth filter on session ' session]); 
        [B,A]   = butter(order,cutoff);              % 6th order lowpass Butterworth filter with maximal cutoff frequency                   
        volumes = filter(B,A,volumes,[],2);        % Filter along columns, row by row (along time)
    elseif strcmp(ftype,'golay')
        % Create a savitsky golay filter, i.e. polynomial fit-smoothing
        disp(['Running Savitsky-Golay filter on session ' session]); 
        b = sgolay(order,windowSize);
        volumes = sgolayfilt(volumes,order,windowSize,[],2); % Filter along the columns, row by row (along time)
    end
    disp('Filtering complete.')
    
    % Add mean back     
    volumes = volumes + (repmat(volume_means,1,reshaped_vol_dims(2)));
    
    % If display is desired, save one voxel for display purposes
    if display; voxel = 1001; post_filter_signal = volumes(voxel,:); end    
    
    % Reshape images back into original format         
    volumes = reshape(volumes,vol_dims);      % Shape back into 4D (you can verify correct operation by not filtering)                              
    
    
    % Save the images. Proceed image-wise, to not mess up different headers
    done = 0; progressbar(['Writing out filtered images from session: ' session]) 
    for I = 1:numel(vol_info)
        % Update file header
        new_vol_info(I) = vol_info(I);
        new_vol_info(I).descrip = [vol_info(I).descrip  ' - ' ftype ' temporal smoothing'];
        new_vol_info(I).private.descrip = new_vol_info(I).descrip;
        %
        new_vol_info(I).fname = [folder,filesep,'f',filenames{I}];
        new_vol_info(I).private.dat.fname = [folder,filesep,'f',filenames{I}];                
        
        % Write out the image
        image = squeeze(volumes(:,:,:,I));
        V(I) = spm_write_vol(new_vol_info(I),image);
        
        done = done+1; progressbar(done/numel(vol_info));
    end     
        
    % If needed, show example voxel pre and post filtering, and events
    if display
        if (F==1)
            filter_figure = figure('position',[0 0 1600 250*numel(folders)]);
        end        
        nrows = numel(folders);        
        row   = F;
        plotnumber = (F*2)-1;
    
        % Prepare event markers if they exist
        timeaxis = 0:TR:(TR * vol_dims(4))-TR; 
        if exist('task_frequency'); task_rate = 1/task_frequency; task_events = 0:task_rate:round(max(timeaxis));  end        
        
        % Create the plots              
        subplot(nrows,2,plotnumber); 
        plot(timeaxis,pre_filter_signal); title(['Session ' session ' before temporal filtering']);
        if exist('task_frequency'); text(task_events,repmat(0.9*max(pre_filter_signal),numel(task_events),1),'X'); end
        xlabel('Seconds in session'); ylabel('Neural activity (a.u.)')
        box off
        %
        subplot(nrows,2,plotnumber+1); 
        plot(timeaxis,post_filter_signal); title(['Session ' session ' after temporal filtering']);
        if exist('task_frequency'); text(task_events,repmat(0.9*max(post_filter_signal),numel(task_events),1),'X'); end
        xlabel('Seconds in session'); ylabel('Neural activity (a.u.)')               
        box off
                
        % Save figure
        saveas(filter_figure,[folder,filesep,'Temporal ' ftype ' filtering.png'],'png')
    end
    
    % Clear large variables
    clear volumes reshaped_volumes dataOut backshaped_images image V new_vol new_vol_info
    disp(['Filtering ' session ' took ' num2str(filtertic) ' seconds.'])
    
end

toc_call = toc(tic_call);
disp(['Temporal filtering completed after ' num2str(toc_call/60) ' minutes.'])
disp('==========================================================')

end
