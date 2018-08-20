function [combined_roi_list] = combine_index_roi_lists(lists)
%% Combine voxel index lists
% Combine voxel index lists simply by appending them.
%
% Input: 
%    - lists: cell string array with Nx1 voxel list vectors
%
% Output:
%    - nx1 vector of combined mask list
%    - dims: dimensions of the corresponding space 
%
% C. Utzerath, 2014-15


%% Combine lists
temp = [];
for i = 1:numel(lists)    
    temp = [temp; lists{i}];
end

combined_roi_list = temp;

end

