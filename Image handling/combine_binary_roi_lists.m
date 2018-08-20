function [combined_roi_list] = combine_binary_roi_lists(lists)
%% Combine binary voxel lists
% Reads in a set of binary roi lists and combines them.
% Result is a binary list of equal length as the input lists, but with
% overlapping ones.
%
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
temp = zeros(numel(lists{1}),1);
for i = 1:numel(lists)    
    temp = temp | lists{i};
end

combined_roi_list = temp;

end

