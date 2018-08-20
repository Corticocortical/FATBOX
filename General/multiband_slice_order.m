function [slice_order] = multiband_slice_order(mb_factor,n_slices,reference)
% Calculate the slice order for a multiband sequence. By default, slice 1
% is the reference slice. Note that correct slice order is also verified
% for when reference slice is 'bottom'; I did not test the other way
% around so proceed with due caution.
%
% Input arguments:
% - mb_factor: the multiband acceleration factor, e.g. 4
% - n_slices: total number of slices.
% - reference: string. 'top' or 'bottom' (default if left empty)
%
% Explanation courtesy of Jose: %
% "If you divide your number of slices by MB factor you get the spacing between simultaneously excited slices
% In the case where this number is 9
%  
% slices 1 10 19 28  are acquired at the same time
% plus "interleaved",
% and slice 1 is the first to be acquired, the slice ordering is:
% 
% [1:9:36 3:9:36, 5:9:36, 7:9:36, 9:9:36, 2:9:36,  4:9:36, 6:9:36, 8:9:36]
%  
% this if the spacing between simultaneously excited slices is odd, otherwise you would start with the even numbers"
%
% Christian Utzerath 2017 (Donders Institute)


%% Parse input
mb = mb_factor;     % Acceleration factor
ns = n_slices;      % N slices
sp = ns/mb;         % Spacing between slices

% Set reference - default is bottom
if strcmp(reference,'bottom')
    rf = 'bottom';
elseif strcmp(reference,'top')
    rf = 'top';
else
    rf = 'bottom';
end
    

%% Calculate
% Determine the slices you walk through without the ones simultaneously
% excitated
if isEven(sp);
    slices = [1 2:2:ns]; 
else
    slices = [1:2:ns];
end

% Now add simultaneously excitated slices
slo = [];
for s = 1:sp
    slo = [slo s:sp:ns];    
end
 
% If you start at the top, flip the order around
switch rf
    case 'top'
        slice_order = fliplr(slo);
    case 'bottom'
        slice_order = slo;
end


end

