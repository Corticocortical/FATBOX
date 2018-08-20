function[] = rsa_plot_RDMs(RDMs,names,cmap,con_names,headline)
% Plots RDMs nicely. Note that RMDs should all be of the same type (either
% correlation or distance), otherwise scaling and color map will be wrong
% since these are all based on the first RDM. Whether or not an RDM shows
% distance or correlations is determined based on the RDM's diagonal.
%
% Names for RDMs, conditions, and figure asl well as plot title can either
% be left open or generated automatically. Function will attempt to use
% brewer maps if no color map is specified. If brewer maps can't be found
% on MATLAB search paths, function defaults to jet and hot as default
% schemes for correlation and distance.
%
% Recommended functions to be on path (helper functions should be on path
% anyways, so should not be a problem)
% - brewermap
% - numSubplot 
% 
% 
% Input arguments:
%       RDMs: Cell array consisting of matrices.
%       Names: Cell string array. If left empty, will use generic name.
%       cmap: string or matrix specifying colormap. If left empty, script
%       will pick based on whether RDM contains distance or correlation.
%       con_names: cell string array with names of conditions. if left
%       empty, will just label them A,B,C,D...
%       headline: string with title for figure
%
% C. Utzerath, 2015

%% Create random RDMs and reset names, cmap (for testing)
% N = 10;
% names = {};
% cmap = '';
% con_names = {}; 
% headline = '';
%
% for i = 1:10;
%     d = rand(N,1); % The diagonal values
%     t = triu(bsxfun(@min,d,d.').*rand(N),1); % The upper trianglar random values
%     RDMs{i} = diag(d)+t+t.'; % Put them together in a symmetric matrix
%     
%     % Set diagonale to 0 (or 1 if you want correlations)
%     for m = 1:N
%        RDMs{i}(m,m) = 0;
%     end    
% end

%% Deal with input
% Check whether RDMs are stored as cell, and if so what format
if iscell(RDMs)   
else
    disp('rsa_plot_RDMs: RDMs must be supplied in a cell. Abort plotting.');
    return
end

% Correlation or distance?
if RDMs{1}(1,1) == 1
    isCorr = 1; isDist = 0;
    style = 'divergent';
else
    isCorr = 0; isDist = 1;
    style = 'sequential';
end

% Check if (sufficient) names for RDMs are supplied or invent some
if ~isempty(names)
    if numel(RDMs) > numel(names)
        for i = [numel(names)+1 : numel(RDMs)]
            names{i} = ['RDM' num2str(i)];
        end
    end
else
    for i = 1:numel(RDMs)
        names{i} = ['RDM ' num2str(i)];
    end
end

% Provide automatic labels for unlabeled conditions
if ~isempty(con_names)
    if size(RDMs{1},1) > numel(con_names)
        for i = [numel(con_names)+1 : size(RDMs{1},1)]
            con_names{i} = ['C' num2str(i)];                
        end
    end
else
    for i = 1:size(RDMs{1},1)
        con_names{i} = ['C' num2str(i)];
    end
end
        

%% Actually plot
% Clumsily determine amount of subplots (neatly if numSubplots is on path)
nplot = numel(RDMs)+1; % reserve extra panel for legend
try
    plotdims = numSubplots(nplot);
catch % if numSubplots isn't on path
    if nplot == 1
        plotdims = [1 1];
    elseif nplot < 5
        plotdims = [2 2];
    elseif nplot < 9
        plotdims = [2 4];
    elseif nplot < 13
        plotdims = [3 4];
    elseif nplot == 15
        plotdims = [3 5];
    elseif nplot == 16
        plotdims = [4 4];
    elseif plotdims < 21
        plotdims = [4 5];
    elseif plotdims < 26
        plotdims = [5 5];
    else
        disp('Could not find numSubplots.m on path and amount of RDMs exceeds hardcoded list of plot arrangements.')
        return
    end
end


% Plot
figure
for iplot = 1:nplot   
    % Nplot refers to last slot in figure, which is legend
    if iplot < nplot      
        % Scale if RDM shows distance (easier comparison between subjects)
        if isDist
            data = rescale_matrix(RDMs{iplot});
        else
            data = RDMs{iplot};
        end
        
        % Plot
        subplot(plotdims(1),plotdims(2),iplot)
        imagesc(data)
        title(names{iplot})
        set(gca,'XTick',[1:numel(con_names)]);
        set(gca,'YTick',[1:numel(con_names)]);
        set(gca,'XTickLabel',con_names)
        set(gca,'YTickLabel',con_names)        
    end
    
    % Last plot is  legend. Note that distance and correlations have
    % different scale and thus different color bars
        if iplot == nplot
            subplot(plotdims(1),plotdims(2),nplot)
            if isDist
                imagesc(linspace(1,0,100)')
            else
                imagesc(linspace(1,-1,100)')
            end            
            
            title('Colormap')
            set(gca,'XTickLabel',{});
            set(gca,'YTick',[1 100]);
            set(gca,'YtickLabel',{'max' 'min'});
        end    
end


% Choose an appropriate color map. Try user-supplied, else decide based on
% measure. Fall back on built-in ones in case something goes wrong.
if ~isempty(cmap)
    try
    colormap(cmap)
    catch
        disp('Custom color map does not work, using my own.')
        if isCorr
            cmap = 'jet' ;
        else
            cmap = 'hot';
        end
    end
else
    switch style
        case 'sequential'
            try
                cmap = brewermap(100,'*BuPu');              
            catch
                disp('Could not find brewer maps, using built-in ones.')
                cmap = 'hot';
            end
        case 'divergent'
            try
                cmap = brewermap(100,'*RdYlBu');
            catch
                disp('Could not find brewer maps, using built-in ones.')
                cmap = 'jet';
            end
    end
    colormap(cmap);
end
    
    
% Add a title
if ~isempty(headline)
    suptitle(headline)   
else
    if isDist
        suptitle('Distance RDMs (rescaled to [0 1])')
    else
        suptitle('Correlation RDMs')
    end
end






end