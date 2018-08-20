function [conjunction] = make_conjunction_mask(list1,list2,save_location)
% Creates a conjunction ROI mask based on two ROI lists. 
%
% Input:
% list 1: nx1 binary vector, with ones for included voxels
% list 2: nx1 binary vector, with ones for included voxels
% save_location: full file path indicating where to store result. If empty,
% no saving.
%
% Output:
% conj: nx1 binary vector with the overlap between both
%
% C. Utzerath, 2015

cm   = zeros(numel(list1),1);
conj = find(list1 == 1 & list2 == 1);
cm(conj) = 1;

conjunction = cm;


% If save, write image at location
if isempty(save_location)
else
    save(save_location,'conjunction')
end


end

