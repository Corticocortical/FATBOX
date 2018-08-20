function [fc,fc_norm] = SG_filter_cutoff(orders,frames,SF,display)
% For a given combination of order and frame, calculates corresponding
% cutoff frequency in Hz. 
%
% Function accepts vectors as inputs for oder and frames, and can hence
% output a matrix giving the approximate cutoff frequency for combinations
% of different polynomial orders and frames.
% 
% Uses formula demonstrated in Schafer 2011: What is a Savitzky-Golay
% filter? IEEE Signal Processing Magazine, 28(4), 111-117.
%
% Note that these equations are only approximate. Below fames < 25 and
% surely frames < 10, the formula becomes less exact.
%
%
% Input arguments:
% - orders: scalar/vector. Polynomial orders to consider.
% - frames: scalar/vector. Window sizes to consider for SG filter.
% - SF: sclar. Sampling frequency (e.g., 1/TR).
% - display: logical. If 1, will show plot.  
%
% Output arguments:
% - fc: frames x orders matrix containing cutoff frequencies in Hertz.
% - fc_norm: frames x orders matrix containing cutoff frequencies in w/pi
%       (i.e., normalized to SF).
%
%% Warn
if numel(find(frames<25))>0
    disp('. Warning: predicted cutoff of SG filter can be imprecise if window size < 25 (see doc)')
end

%% Calculate
fc = zeros(numel(frames),numel(orders));
fc_norm = zeros(numel(frames),numel(orders));
for o = 1:numel(orders)
    for f = 1:numel(frames)
        % Constant term (see equation 12 and corresonding text in Schafer 2011)
        if frames(f) < 10
            constant_term = 2;
        elseif frames(f) >= 10 && frames(f) < 25
            constant_term = 2;
        elseif frames(f) >= 25
            constant_term = 4.6;
        end
        
    fc_norm(f,o) = (orders(o) + 1) / ((3.2*frames(f) - constant_term));
    fc(f,o) = fc_norm(f,o)*SF;                                
    end
end

%% Display
if display
    % Get dimensionality of data to see if matrix or vector (ismatrix /
    % isvector don't do this nicely)
    result = figure('position',[0 0 800 400]);
    dims = size(fc);
    if dims(1) > 1 && dims(2) > 1
       imagesc(fc)
       title('Predicted SG filter cutoffs (Hz) as function of polyn. order/window size','fontsize',12)
       xlabel('Polynomial order','fontsize',12); ylabel('Window size','fontsize',12);       
       set(gca,'xtick',1:numel(orders)); set(gca,'xticklabel',orders)
       set(gca,'ytick',1:numel(frames)); set(gca,'yticklabel',frames)
       colorbar
       box off
    elseif dims(1) == 1 && dims(2) == 1
        disp('No plot being made because only one fc computed.')
    else
       plot(fc)  
       xlabel('Parameter varied (polynomial order or window size)','fontsize',12); ylabel('Cutoff frequency (Hz)','fontsize',12);
       title('Predicted SG filter cutoff frequencies as function of polynomial order/window size','fontsize',12)
       box off
    end
    
  


end

end