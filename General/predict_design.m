function [D] = predict_design(onsets,names,duration_per_con)
% Shows how the design matrix of a given run might look like 
% (approximately). Useful for checking whether your model works. 
%
% Input argument: 
% onsets: cell matrix array onsets n scans, as SPM would have them
% names: cell string array with condition names
% duration_per_con: vector with the duration of of each condition
%
% Christian Utzerath 2014-15 (Donders Institute)


durs = duration_per_con;
ncon = numel(names);
ind = {};
D   = [];
for c = 1:ncon
    for i = 1:durs(ncon)+1
        ind{c} = onsets{c} + i-1;   
        D(ceil(ind{c}),c) = 1;
    end
end




end

