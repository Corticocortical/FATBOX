function [ output_args ] = qc_plot_subject_motion(folders,limit,save)
% This function shows the subject's motion during successive runs of an
% experiment. It will do so by reading out a parameter file (rp*.txt) that
% is usually generated during realignment with programs such a SPM. 
%
% Input arguments:
% - Folders: cell with strings leading to the folder(s) in which the images
%            of each session are stored along with the rp file(s).
% - Limit:   scalar. Define a limit in mm that, when exceeded, triggers a
%            warning.
% - Save:    Binary scalar. 1: save an image of the subject's head motion.
%            Will be saved in a folder at the same level as the sessions.
%            Note that all <folders> are assumed to have the same level in
%            the folder hierarchy.
%
% NEEDS: 
% - realignment parameters in each specified folder. 
%
% Christian Utzerath 2015 (Donders Institute)


%% Do admin things
% Find and/or create a quality check folder at the level of <folders>
return_dir = pwd;
cd(folders{1})
cd ..
toplevel = pwd;
qc_dir = [toplevel,filesep,'Quality'];
if ~exist(qc_dir); mkdir(qc_dir); end
cd(return_dir)

% Read in motion
RPs = {}; 
for f = 1:numel(folders)
    % Find and read file
    cd(folders{f})    
    clear d
    d = dir('rp*.txt');    
    if numel(d) > 1
        disp('Warning: multiple realignment parameter files found. Abort.')
        return
    else
        file = d.name;
        fid = fopen([folders{f},filesep,file],'r');
        RPs{f} = fscanf(fid,' %e %e %e %e %e %e',[6,inf])';
        fclose(fid);        
    end
end
cd(return_dir)

% %% Plot (translation + rotation)
% % Determine number of plots and canvases needed (1 row per RP file)
% n         = numel(RPs);
% n_per_fig = 4;
% n_fig     = ceil(n/n_per_fig); 
% 
% % Open first figure, start counting plots
% figure('position',[0 0 900 600])
% 
% % Make index numbers for plots to divide on left and right half
% lefts  = 1:2:100;
% rights = 2:2:100; 
% if n > max(lefts)
%     disp('Too many plots too draw!')
%     return
% end
% 
% % Start drawing sessions. If 4 are drawn, open new figure.
% ss = 1; % keeps track of sessinons drawn on this particular figure
% for s = 1:numel(RPs)
%     % Session name
%     session = folders{s};
%     dashpos = max(strfind(session,filesep));
%     sname   = session(dashpos+1:end);
%     
%     % Data
%     motion   = RPs{s}(:,1:3);
%     rotation = RPs{s}(:,4:6); 
%     
%     % Translation plot    
%     subplot(4,2,lefts(ss));
%     plot(motion);
%     title(['Translation during ' sname]);
%     legend({'x' 'y' 'z'});
%     xlabel('Scans'); ylabel('mm');
%     
%     
%     % See if there need to be triggers/warnings added
%     clear triggers
%     triggers = find(motion > limit);
%     if ~isempty(triggers);
%         title(['EXCESSIVE MOTION during ' sname]);
%     else
%          title(['Translation during ' sname]);
%     end
%     
%     % Rotation 
%     subplot(4,2,rights(ss));
%     plot(rotation);
%     title(['Rotation during ' sname]);
%     legend({'a' 'b' 'c'});
%     xlabel('Scans'); ylabel('degrees');
%     
%     % Update counter for plot row
%     ss = ss +1;
%     
%     % Open new figure once four sessions are drawn
%     if rem(s,4) == 0        
%         % Save the old figure first
%         saveas(gcf,[qc_dir,filesep,'Subject motion until session ' num2str(s),'.png']);
%         
%         figure('position',[0 0 900 600])
%         ss = 1;
%     end
%     
%     
% end
% 
%      %Save the last figure as well
%     saveas(gcf,[qc_dir,filesep,'Subject motion until session ' num2str(s),'.png']);

%% Alternative plot (translation only)
% Open  figure
figure('position',[0 0 900 600])

% Start drawing sessions.
for s = 1:numel(RPs)
    % Session name
    session = folders{s};
    dashpos = max(strfind(session,filesep));
    sname   = session(dashpos+1:end);
    
    % Data
    motion   = RPs{s}(:,1:3);
    rotation = RPs{s}(:,4:6);
    
    
    % Plot only translation
    subplot(3,2,s)
    plot(motion)
    title(['Translation during ' sname]);
    legend({'x' 'y' 'z'});
    xlabel('Scans'); ylabel('mm');
    
    % See if there need to be triggers/warnings added
    clear triggers
    triggers = find(motion > limit);
    if ~isempty(triggers);
        title(['EXCESSIVE MOTION during ' sname]);
    else
        title(['Translation during ' sname]);
    end    
end


% Save figure
saveas(gcf,[qc_dir,filesep,'Subject head motion.png']);




        

end

