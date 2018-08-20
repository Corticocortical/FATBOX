function[] = rsa_simple_distplot(D,names,p_title,colormapping,dotcolor,dotsize,axislabels,fontsize_chart,fontsize_dotlabels)
%% Simple distances plot
% Uses classical multidimensional scaling to visualize distance between
% different conditions.
%
% Uses an undocumented feature: line smoothing. Feature has the bug that
% smoothed lines will obstruct dots and annotation. Recommended workaround
% is to have their location specified in 3D, and increase the z-coordinate
% by 1. Then they are above the smoothed lines in z-dimension and are shown
% again.
%
%
% Input argument: 
%   - D: square, symmetric distance matrix
%   - names: cell string array with names of conditions (can be left empty)
%   - p_title: string with title for plots, can be left empty.
%   - colormapping: custom color map. If left empty ([]), generates its own.
%   - dotcolor: vector that assigns each condition/dot a color from the 
%     color map. If left empty ([]), all conditions/dots get same color.
%   - dotsize: vector that specifies size of each dot. If left empty,
%     defaults to 40.
%   - axislabels: the label that is put on x and y axis. Empty if left
%     empty.
%   - fontsize_chart: font size for title and labels. If empty, defaults.
%   - fontsize_dotlabels: font size for dot labels. If empty, defaults.
%
%
% C. Utzerath, 2014-15

%% Deal with input arguments
% Check if D has correct format
if size(D,1) == size(D,2)
else
    disp('rsa_simple_distplot: matrix has improper format!')
    return
end

% Improvise names if necesary
if isempty(names); 
    names = {};
    for i = 1:size(D,2)
        names{i} = num2str(i);
    end
end

D = D;
names = names;


% Check title
if strcmp(p_title,'')
    p_title = 'Distances between conditions projected onto 2D space';
end

% Check if color argument is supplied for dots, and if not set default
if isempty(dotcolor)
    colormapping = [137 204 240];   
    colormapping = colormapping./255;
    color_correspondence = repmat(ones,numel(names),1);
else
    colormapping = colormapping;
    color_correspondence = dotcolor;
end

% Check if dot size supplied
if isempty(dotsize)
    dot_size = repmat(40,numel(names),1);
else
    dot_size = dotsize;
end


% Set font sizes
if isempty(fontsize_chart)
    fontsize_chart = 16;
end
if isempty(fontsize_dotlabels)
    fontsize_dotlabels = ceil(dot_size/8);
end



%% Program
try
    % Let matlab compute classical multidimensional scaling
    [Y,eigvals] = cmdscale(D);
    
    % Look at eigenvalues
    format short g;
    [eigvals eigvals/max(abs(eigvals))];
    
    % Do magic
    Dtriu = D(find(tril(ones(size(D,2)),-1)))'; % triangular form
    maxrelerr = max(abs(Dtriu-pdist(Y(:,1:2))))./max(Dtriu);
    
    
    % Plot the dots for each conditions projected on 2D space
    % Create lines between data points (credit to gnvoice / Stackoverflow)
    grey = [0.9,0.9,0.9];                        % RGB value for grey.
    N = length(Y);                               %# Number of points
    x = Y(:,1);                                  %# Position of dots on x axis
    y = Y(:,2);                                  %# Position of dots on y axis
    [r,c,v] = find(hankel(2:N));                 %# Create unique combinations of indices
    index = [v c].';                             %'# Reshape the indices
    line(x(index),y(index),'Color',grey,'linesmoothing','on');        %# Plot the lines       
        
    % Prepare 3D scatter plot for dots (3D necessary due to line smoothing)    
    hold on
    z = ones(length(x),1);                                              % Fake third dimension
    distplot = scatter3(x,y,z,dot_size,color_correspondence,'filled');  % Add data points, filled
    set(distplot,'MarkerEdgeColor','none')                                 % Give points edge
    colormap(colormapping)                                              % Switch to your colormap so color_correspondence matches to anything
    hold off
    
    % Add labels a bit to the side of dots (play around with these numbers)
    offset = range(Y(:,1) * 0.01 * dot_size/40);    
    h=text(Y(:,1)+offset, Y(:,2), ones(length(Y),1),names)      % Put condition names next to the points. Mind you: line smoothing requires text to be placed on a z-dimension!   
    set(h,'FontSize',fontsize_dotlabels);
    
    % Style the axes and title         
    xlim([min(x)-offset*5 max(x)+offset*5]); ylim([min(y)-offset*5 max(y)+offset*5]);    
    set(gca,'xtick',[]); xlabel(axislabels,'FontSize',fontsize_chart)
    set(gca,'ytick',[]); ylabel(axislabels,'FontSize',fontsize_chart)
    title(p_title,'FontSize',fontsize_chart)       
    

    
    
catch
    disp('Multidimensional scaling and subsequent plotting were not possible for some reason,')
    disp('Possibly, distance matrix has too many negative eigenvalues?')
end



end





