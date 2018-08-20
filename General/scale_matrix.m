function [scaled_M] = scale_matrix(M,dim,mode)
% Scales a matrix, such as a data matrix. 
%
% Input arguments: 
% - M: matrix to scale. 2D numerical.
% - dim: numerical, 1, 2, or 2, specifying dimensions.
%                 * 1: scale across rows, giving every column the same mean.
%                 * 2: scale across columns, giving every row the same mean.
%                 * 3: scale all matrix elements togehter.
% - mode: string. scaling mode. Options are: 
%                 * '01'    : sets the maximum value to 1 and the minimum
%                             to 0 (see matlab mat2gray) (default)
%                 * 'z'     : uses Matlab's zscore function to zscore
%                 * 'sum'   : scales every element to the magnitude of
%                             sum of elements: element./abs(sum(elements))
%                             Note: this does not guarantee same mean for
%                             rows and columns!
%                 * 'mean'  : divides elements by the mean (scaling prop.
%                             to mean)
%                 * 'none'  : do nothing, output = input 
%
% This function is not optimized and works particularly well on small matrices.
%
% Christian Utzerath 2015 (Donders Institute)


%% Check input arguments
if ~ismatrix(M)
    disp('Error in scale_matrix: input provided is not a matrix,')
    return
end

if isempty(dim)
    dim = 2;
end

if isempty(mode)
    mode = '01';
end

%% Do nothing
if strcmp(mode,'none')
    rM = M;
end

%% Greyscale scaling
if strcmp(mode,'01')
    rM = [];
    switch dim
        case 1
            for c = 1:size(M,2)
                rM(:,c) = mat2gray(M(:,c));
            end
        case 2
            for r = 1:size(M,1)
                rM(r,:) = mat2gray(M(r,:));
            end
        case 3
            rM = mat2gray(M);
    end
end

%% Z scaling
if strcmp(mode,'z')
    rM = [];
    switch dim
        case 1
            rM = zscore(M,0,1);
        case 2
            rM = zscore(M,0,2);
        case 3
            m  = reshape(M,numel(M),1);
            rm = zscore(m);
            rM = reshape(rm,size(M,1),size(M,2));
    end
end

%% Sum scaling
if strcmp(mode,'sum')
    rM = [];
    switch dim
        case 1
            for c = 1:size(M,2)
                %denominator  = abs(sum(M(:,c)));
                denominator = sum(M(:,c));
                rM(:,c) = M(:,c)./denominator;
            end
        case 2
            for r = 1:size(M,1)
                %denominator = abs(sum(M(r,:)));
                denominator = sum(M(r,:));
                rM(r,:) = M(r,:)./denominator;                
            end            
        case 3
            %denominator = abs(sum(sum(M)));
            denominator = sum(sum(M));
            rM = M./denominator;
    end
end

%% Mean scaling
if strcmp(mode,'mean')
    rM = [];
    switch dim
        case 1
            for c = 1:size(M,2)
                denominator  = sum(M(:,c)) / numel(M(:,c));
                rM(:,c) = M(:,c)./denominator;
            end
        case 2
            for r = 1:size(M,1)
                denominator = sum(M(r,:)) / numel(M(r,:));
                rM(r,:) = M(r,:)./denominator;
            end            
        case 3
            denominator = mean(M);
            rM = M./denominator;
    end
end

%% Output
scaled_M = rM;