%% Setup our little filtering experiment
% This file is an attempt to simulate, and then build with real data, what
% happens when we high-pass filter fMRI data. Note the following points of
% thought that need dealing with: 
%
% 1. Bugs
% In the plots, some of the time axis might be off (confusion due to
% different sampling frequencies in SPM's HRF and our fMRI data).
%

% We are using SPM, so put it on path
add_SPM;

% Parameters
TR             = 0.68;     
task_rate      = 14;          % Rate of stimulus presentation (s) (trial length)
task_frequency = 1/task_rate; % Frequency of task is about 1/13.6 seconds 
n_trials        = 10;         % Simulate N trials
trial_dur       = 14;         % How long should an event last
card_frequency = 0.5;         % Possible frequency of heart rate

% Select filter ('moving average' or 'lowpass')
filtertype = 'moving average';

%% Create simulated physiological signal
[hrf,p,hrf_sampling_rate,blood,bloody_hrf,noisy_bloody_hrf] = simulate_noisy_HRF(TR,card_frequency,0.1,0.1,0);

%% Simulate a voxel's time course in the experiment
total_dur = n_trials * task_rate;   % approx. dur. of simulated trials (excluding what is added when convoluting signals)
n_repeats = total_dur/p(7);         % p(7) is duration of hrf in s

% Place events in a vector
onsets      = 1:task_rate:total_dur;                    % Trial onsets (s)
event_vec   = zeros(1,total_dur*hrf_sampling_rate);     % Time samples in simulation. Note that sampling rate is specified in sampling rate of spm's simulated HRF, not TR!
timeaxis    = linspace(0,total_dur,numel(event_vec));   % Time points in simulation
for tr = 1:n_trials
   % Assign every trial onset to the closest sampling point for convenience
   onset = onsets(tr);
   difference = timeaxis - onset;         % eg, if a sampling point fell together with an onset, difference would be 0
   nearest_sample = min(abs(difference)); % find nearest sample
   idx_sample     = find(abs(difference) == nearest_sample); % corresponding index         
   
   % Put a stick in the events vector there to model stimulus onset
   % Also set subsequent samples to 1 to simulate 
   length_window = ceil(trial_dur / hrf_sampling_rate);
   event_vec(idx_sample:idx_sample+length_window-1) = 1;
end

%% Create and show signal
% Convolute HRF with event series. Create time vector for all samples
signal   = conv(event_vec,noisy_bloody_hrf);
timeaxis = 0:TR:TR*(numel(signal)-1);

% Add a bit more noise so the signal isn't the same everywhere
signal = signal + 0.35*rand(1,numel(signal));

% Signal without any noise, i.e. ground truth
ground_truth = conv(event_vec,hrf);

%% Design filter
if strcmp(filtertype,'lowpass')
    % Using a low pass
    cutoff = 1/TR/2;                    % max. cutoff frequency is half the sampling rate (cf. Nyquist). That's just above our task rate, so, okay!
    [B,A] = butter(6,cutoff);           % 6th order lowpass Butterworth filter with maximal cutoff frequency
    %freqz(B,A)                         % Show impulse response of filter
    dataIn = signal;
    dataOut = filter(B,A,dataIn);
elseif strcmp(filtertype,'moving average')
    % Using a moving window
    windowSize = 5;                                     % Three images at a time (3*0.68=2.04s)
    b = (1/windowSize)*ones(1,windowSize); a = 1;       % See formula/example doc filter
    dataIn = signal;
    dataOut = filter(b,a,signal);        % Filter along columns, row by row (along time)
end
%% Show overview
figure; 
nrows = 6;
time = linspace(0,p(7),numel(hrf));
subplot(nrows,2,1)
plot(time,hrf); 
title(['Hemodynamic response function (TR =  ' num2str(TR) ' s)'] );
xlabel('Time (s)'); ylabel('Voxel signal (a.u.)')
%
subplot(nrows,2,2)
plot(time,blood); 
title('Blood signal contribution (raw)')
xlabel('Amplitude'); ylabel('Time (s)')
%
subplot(nrows,2,3)
plot(time,bloody_hrf); 
title('HRF modulated by noised blood fluctuations')
xlabel('Time (s)'); ylabel('Voxel signal (a.u.)');
%
subplot(nrows,2,4)
plot(time,noisy_bloody_hrf); 
title('HRF modulated by blood and random noise')
xlabel('Time (s)'); ylabel('Voxel signal (a.u.)');
%
subplot(nrows,2,5:6)
plot(timeaxis,ground_truth)
xlabel('Seconds during experiment'); ylabel('Voxel signal (a.u.)');
title(['Ground truth: neural response.'])
text(onsets,ones(numel(onsets),1),'X')
%
subplot(nrows,2,7:8)
plot(timeaxis,signal)
xlabel('Seconds during experiment'); ylabel('Voxel signal (a.u.)');
title(['Voxel activity over time (TR = ' num2str(TR) '), HR = ' num2str(card_frequency) ', ' num2str(n_trials) ' trials.'])
text(onsets,ones(numel(onsets),1),'X')
%
subplot(nrows,2,9:10)
plot(timeaxis,dataOut)
xlabel('Seconds during experiment'); ylabel('Voxel signal (a.u.)');
text(onsets,ones(numel(onsets),1),'X')
title(['Voxel activity after butter filter (cutoff = ' num2str(cutoff) 'Hz)'])
%
subplot(nrows,2,11:12)
c1 = corrcoef(ground_truth,dataIn); c1 = c1(1,2);
c2 = corrcoef(ground_truth,dataOut); c2 = c2(1,2);
bar([c1 c2]); minima = min([c1 c2]); maxima = max([c1 c2]);
ylim([(minima-0.05*minima)   (maxima+0.05*maxima)]);
set(gca,'xticklabel',{'Measured signal' 'Cleaned signal'})
xlabel('Signal type')
title('Corelation between filtered signal and ground truth')
%
suptitle(['Voxel simulation for ' filtertype ' filter'])
