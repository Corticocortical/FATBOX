function [amplitude_spectrum] = qc_amplitude_spectrum_volumes(volumes,plot_options)
% This function looks for spiking slices (elevated mean activity) across a
% series of volumes. No report will be saved.
%
% Input: 
% - volumes: cell with strings, each being a path to a volume.
% - plot_options: 1xn cell array with flags for optional plotting. Elements:
%               {1}: binary. 1: show a plot of the amplitude spectrum
%               {2}: binary. 1: make a new figure (if empty defaults to 0
%                    so you can plot into existing figures).
%               {3}: colormap. either a string for a default colormap or a 
%                    nx3 RGB color matrix.
%               {4}: show colorbar. 1: show colorbar.
%               If any flag is empty, function uses defaults.
%
% Output: 
% - amplitude_spectrum: a volume x slice matrix with activity values
%
% Christian Utzerath 2015 (Donders Institute)


%% Parse options
% Check whether all arguments are present or give defaults
if isempty(plot_options{1})
    plot_spectrum = 1;
else
    plot_spectrum = plot_options{1};
end

if isempty(plot_options{2})
    new_fig       = 0;
else
    new_fig       = plot_options{2};
end

if isempty(plot_options{3})
    map           = 'jet';
else
    map           = plot_options{3};
end

if isempty(plot_options{4})
    show_map      = 1;
else
    show_map      = plot_options{4};
end


% Determine colors - either load a map, or use the existing one
if ischar(map)
    map = map;
elseif ismatrix([map]) && size(map,2)==3;
    map = map;
else
    disp('Error: invalid color map specified.')    
end

%% Process
% Get amplitude spectrum
amplitude_spectrum = [];
for v = 1:numel(volumes)
   amplitude_spectrum(v,:) = qc_amplitudes_slices(volumes{v});     
   
   % Show you're not dead 
   if numel(volumes) > 100 && ~rem(v,50)
        disp(['Computing amplitude spectra: ' num2str(v/numel(volumes)*100) '% ' ])
   end
end

% Show spectrum
if plot_spectrum
    
    if new_fig
        figure;
    end
    
    imagesc(amplitude_spectrum);
    xlabel('Slices')
    
    
    % Y axis 
    ylabel('Volumes')
    %ticks    = ceil(linspace(1,numel(volumes),10))
    %set(gca,'YTick',1:numel(volumes))
    
    % Coloring
    colormap(map)    
    if show_map
        colorbar;
    end
    
    title('Spectrum of average amplitude deviations from slice mean in volume series')
    box off
end

end

