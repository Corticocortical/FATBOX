function [window_size,effective_cutoff] = moving_average_window_size_lowpass(TR,task_rate,desired_cutoff,display)
% Given an experiment's TR and the rate of the task, this function
% determines what the most fitting window size is for a moving average
% low-pass filter.
%
% First it predicts how different window sizes will result in different
% frequencies being able to being modelled (i.e. estimate the Nyquist
% frequencies corresponding to different window sizes).
%
% Then it choses the window size (based on the Nyquist frequencies) that
% falls closest to the desired cutoff.
%
% Please note that this is not a cutoff frequency in the sense that might
% be used by electrical engineers - the function does not determine the
% dropoff in dB at particular frequencies, but just estimates what the
% highest frequencies are that you might want to reasonably include.
%
% Input arguments
% - TR: scalar. TR of your experiment (s), e.g. 0.68.
% - task_rate: scalar. Only useful for blocked design, leave empty
%              otherwise. Will report in case task rate is falling above
%              the cutoff.
% - desired_cutoff: scalar. Frequency in Hertz above which contributions
%              should be attenuated from the data.
% - display: logical. Plot Nyquist frequencies as function of window size?
%
%% Do
% Get Nyquist frequencies for different window sizes
for i = 1:20
    TRs(i) = i*TR;
end
Nyquists = convert_TR_to_Nyquist(TRs);

% Pick window size that most closely above desired cutoff
differ = Nyquists - desired_cutoff;
still_above_cutoff   = differ > 0;
least_positive_difference = min(differ(still_above_cutoff==1));
element = find(differ == least_positive_difference);

% Output
window_size = element;
effective_cutoff = Nyquists(element);
disp(['The chosen cutoff frequency is ' num2str(desired_cutoff) 'Hz; the effective cutoff frequency is ' num2str(effective_cutoff) 'Hz.'])

% Verify task rate still works
if effective_cutoff <= 1/task_rate
    disp('WARNING: The effective cutoff is lower than the task frequency; task-related activation will be removed!')
end

% Display
if display
    figure
    plot(TRs,Nyquists)
    xlabel('Window size in TR'); ylabel('Corresponding Nyquist frequencies')
    title('Finding best window size for moving average lowpass for desired cutoff frequency') 
    text(window_size,effective_cutoff,'<- closest window size')
    box off
end


end

