function [Nyquist] = convert_TR_to_Nyquist(TR)
% Simply outputs the Nyquist frequency that you can estimate with your
% given TR. 
%
% Input: 
% - TR: scalar or vector. The TR in the expeirment (in seconds). If TR is a
%       vector, function will be performed on every element in TR.
%
% Output:
% - Nyquist: maximal frequency that you can model at that TR, reasoning
%            that any periodic signal needs to be sampled twice to be 
%            estimated. If TR is a vector, so will Nyquist be.
%
%
%% Do
if numel(TR)==1
    SF = 1/TR;
    Nyquist = SF/2;
else    
    for i = 1:numel(TR)
        SF(i) = 1/TR(i);
        Nyquist(i) = SF(i)/2;
    end
end

end

