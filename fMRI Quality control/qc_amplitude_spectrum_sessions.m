function [amplitude_spectra] = qc_amplitude_spectrum_sessions(folders,prefix,plot_options)
% Uses qc_amplitude_spectrum_volumes to plot the amplitudes of all volumes
% across several sessions. Useful to discover slices or runs that have
% unusually high intensities. 
%
% Will save a report in a quality check folder, which is at the same level
% as the session folders.
%
% Input: 
% - folders: cell string. Folders to scan for volumes
% - prefix: SPM-compatible file filter, e.g. '^arf.*\.nii$'
% - plot_options: 1x3 cell with flags for plotting. Elements/flags are: 
%           {1}: plotting. binary. 1: make (and save) plot.
%           {2}: colormap. string for matlab default or nx3 RGB matrix.
%           {3}: colorbar. binary. 1: show colorbar on plots.
%           If options are empty, default values will be used.
%
%
% Output: 
% - amplitude_spectra: cell array with, per session, amplitude spectrum as 
%                      nvol x nslice array.
%
% NEEDS:
% - SPM to read volumes. 
%
% Christian Utzerath 2015 (Donders Institute)



%% Check and if necessary create quality control folder
return_dir = pwd;
cd(folders{1})
cd ..
toplevel = pwd;
qc_dir = [toplevel,filesep,'Quality'];
if ~exist(qc_dir); mkdir(qc_dir); end
cd(return_dir)

%% Parse input arguments
% Check whether all arguments are present or give defaults
if isempty(plot_options{1})
    plot_spectrum = 1;
else
    plot_spectrum = plot_options{1};
end

if isempty(plot_options{2})
    map           = 'jet';
else
    map           = plot_options{2};
end

if isempty(plot_options{3})
    show_map      = 1;
else
    show_map      = plot_options{3};
end


% Determine colors - either load a map, or use the existing one
if ischar(map)
    map = map;
elseif ismatrix([map]) && size(map,2)==3;
    map = map;
else
    disp('Error: invalid color map specified.')    
end

%% For each session, find volumes and produce amplitude spectra
% Determine number of plots and canvases needed (1 row per RP file)
n         = numel(folders);
n_per_fig = 4;
n_fig     = ceil(n/n_per_fig); 
ss        = 1; % counter
figure('position',[0 0 900 600])

spectra = {};
for f = 1:numel(folders)
    % Identify session
    session  = folders{f};
    fileseps = strfind(session,filesep);
    sessname = session(fileseps(end)+1:end);    
    disp(['Amplitude spectrum for session: ' sessname ])
    
    % Get amplitude for this session
    [volumes] = cellstr(spm_select('FPListRec',session,prefix));
    subplot(4,1,ss);
    spectra{f}  = qc_amplitude_spectrum_volumes(volumes,{plot_spectrum,0,map,show_map});
    title(['Amplitude spectrum for session: ' sessname ])
    
    % Update counter for plot row
    ss = ss +1;
    
    % Open new figure once four sessions are drawn, reset counter
    if rem(f,n_per_fig) == 0        
        % Save the old figure first
        saveas(gcf,[qc_dir,filesep,'Session amplitude spectra until session ' num2str(f),'.png']);
        
        figure('position',[0 0 900 600])
        ss = 1;
    end
end
% Save the last figure as well
saveas(gcf,[qc_dir,filesep,'Session amplitude spectra until session ' num2str(f),'.png']);


amplitude_spectra = spectra;
end

