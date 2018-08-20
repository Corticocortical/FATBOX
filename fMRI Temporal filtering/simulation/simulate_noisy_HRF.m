function [hrf,p,hrf_sampling_rate,blood,bloody_hrf,noisy_bloody_hrf] = simulate_noisy_HRF(TR,card_frequency,noise_card,noise_HF,plotme)
% Simulates a noisy hemodynamic response function. Requires SPM to be on
% path. 
%
% Input: 
% - TR: numeric. the TR of the experiment in seconds.
% - card_frequency: numeric. frequency of cardiac artifacts, i.e. heart frequency
% - noise_card: noise ampl. applied to cardiac artifact (typically < .05)
% - noise_HF: noise ampl. applied to high frequency noise (typically < .05)
% - plotme: logical. 0: show no plot. 1: show plot.
%
% Output: 
% hrf,p - see spm_hrf
% hrf_sampling_rate: sampling rate of HRF
% bloody_hrf: hrf with blood fluctuations added
% bloody_noisy_hrf: hrf with blood fluctuations and added random (HF) noise

%%
% Create a BOLD signal
[hrf,p]           = spm_hrf(TR);
hrf_sampling_rate = numel(hrf) / p(7); 
blood   = sin(card_frequency);

% Create heart rate fluctuations at same sampling rate as SPM's HRF
time = linspace(0,p(7),numel(hrf));
blood = sin(time);

% Apply heart rate fluctuations to HRF, with bit of added randomness
noise_vec  = noise_card*rand(1,numel(hrf));       % Random amplitude for heart noise
bloody_hrf = hrf +  (blood.*noise_vec)';          % Modulate HRF by heart noise

% Add high frequency noise (i.e., just random numbers)
hf_noise = noise_HF*rand(1,numel(hrf));
noisy_bloody_hrf = bloody_hrf + hf_noise';

%% Visualize
if plotme
figure; 
subplot(2,2,1)
plot(time,hrf); 
title(['Hemodynamic response function (TR =  ' num2str(TR) ' s)'] );
xlabel('Time (s)'); ylabel('Voxel signal (a.u.)')
%
subplot(2,2,2)
plot(time,blood); 
title('Blood fluctuations')
xlabel('Amplitude'); ylabel('Time (s)')
%
subplot(2,2,3)
plot(time,bloody_hrf); 
title('HRF modulated by noisy blood fluctuations')
xlabel('Time (s)'); ylabel('Voxel signal (a.u.)');
%
subplot(2,2,4)
plot(time,noisy_bloody_hrf); 
title('HRF modulated by blood and random noise')
xlabel('Time (s)'); ylabel('Voxel signal (a.u.)');
%
suptitle('Simulated voxel')
end



end

