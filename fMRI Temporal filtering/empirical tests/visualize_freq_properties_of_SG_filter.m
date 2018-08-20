% Test effect of filter on signal to figure out whether omega/pi is indeed
% the sampling frequency (i.e. confirm the frequency characteristics of the
% SG filter as outlined by Schaefer)
%
%
% From the Schaefer paper it was a bit unclear whether their plots
% normalize the frequency to Nyquist,or sampling frequemcy, or something
% else. 
%
% They normalize to omega/pi (i.e., omega/pi [0:1])).
%
% After looking up the angular frequency omega, I found it to be 2*pi*f:
% (1) w = 2*pi*f;
% 
% This means that the endpoint of Schaefer's scale might have been :
%
% (2) w/pi = 2*pi*f/pi
%
% By substituting w/pi with 1 we should be able to see what their endpoint
% is:
%
% (3) w/pi = 1 = 2*pi*f/pi
%            1 = 2*f
%
% Which is interesting, because 2*f describes the relationship between
% Nyquist (N) and sampling frequency (SF):
%
% (4) N = SF/2
%
% However, f might also be the sampling rate. Or any other frequency,
% really. But SF and N are the most obious ones. Which one is it? Let's
% find out: 
%   if f == N   --> 1=2*FS means 1 = 2*N = SF
%   if f == SF  --> 1=2*N  means 1 = 2*SF = N   

%% Make a random signal composed of several sine waves
T = 0.68;                                   % 1 sample every xx seconds
SF = 1/T;                                   % sampling frequency
time_axis = 0 : T : T*200;                  % time axis
N  = SF/2;                                  % Nyquist frequency

% Create frequency components
comps = [];
for f = [1:1000]
    frequency(f) = f*0.001;  % ditto 
    %comps(f,:) = sin(2*pi*frequency(f)*time_axis)*2;
    comps(f,:) = rand(100,numel(T))*2;
end

% Add components to complex signal
signal = sum(comps,1);
wn = (rand(numel(time_axis),1)-0.5)*2;
signal = wn;
signal=sinc(time_axis);
plot(time_axis,signal)
xlabel('Seconds'); title('Signal: sinc function')

% Make spectogram of input signal
[freqs,power] = make_power_density_spectrum(signal,SF,T,1);
title('Box-ish power spectrum of input sinc function','Fontsize',12)

%% Filter
% At order 2, farmes 5, signal cutoff (-3db) should be at 0.197*SF, meaning
% the signal should halve there
order  = 2;
frames = 5;
output = sgolayfilt(signal,order,frames);
%
[freqs2,power2] = make_power_density_spectrum(output,SF,T,1);
%

% Where would we expect the cutoff to be, depending on Schaefer's notation?
expected_cutoff = 0.197*SF;
alternative_cutoff = 0.197*N;

% Annotate: mark important frequencies
y1 = get(gca,'ylim');
hold on 
plot([expected_cutoff,expected_cutoff],y1)
plot([alternative_cutoff,alternative_cutoff],y1)
plot([N,N],y1)
hold off
text(expected_cutoff,y1(2)*0.8,'<- cutoff if  w/pi = SF')
text(alternative_cutoff,y1(2)*0.5,'<- cutoff expected if w/pi = Nyquist')
text(N,y1(2)*0.3,'<- Nyquist')
title('Power spectrum after signal was filtered at 0.197*\omega/\pi','fontsize',12)


%% examine a model fitlter
b = sgolay(order,frames);       % B coefficients for our filter
a  = 1;                         % 'The denominator of FIR filters is, by definition, equal to 1' - see help filter
golay = dfilt.df1(b,a);         % design butterworth digital filter object
%
%[b,a] = butter(8,0.25);
fvtool(golay)


%% Calculate cutoff frequency for SG filter  if frames > 25 and order < frames
frames = 7;
order  = 2;
fc = (order+1)/((3.2*frames)-2); % see schaefer paper, for low number of frames
ordinary_fc = fc/SF;


0.165/SF
