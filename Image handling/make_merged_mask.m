function [merged_mask] = make_merged_mask(masks)
% Merges a series of masks into one. You can supply either vectorized 
% binary masks (mask lists), or mask images. 
% 
% Input arguments:
% - masks: cell array. Contains either strings (file paths to mask imagaes) 
%          or nvox x 1 binary vectors (vectorized binary masks / list).
%          Must be consistent for all cells.
%
% Output argument:
% - merged_mask: nvox x 1 binary vector 
%
% C. Utzerath, 2015

%% Deal with input
isImage = [];
if iscellstr(masks)
    isImage = 1;
else
    isImage = 0;    
end

N = numel(masks);
n = numel(masks{1});

%% Program
temp = zeros(n,1);
for i = 1:N
    if isImage==1 
        [cache,dims] = make_list_from_image(masks{i});   % if masks are images, read in first ...            
        if i == 1
            n = numel(cache); temp = zeros(n,1);      % correct the size of temp array when taking in first image
        end
        
        else 
        cache = masks{i};                              % if list, simply take content...
        
    end
    temp = temp | cache;                               % ... and merge.
end
merged_mask = temp;



end

