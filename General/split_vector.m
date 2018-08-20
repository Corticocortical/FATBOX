function [v1, v2] = split_vector (vec)
% Splits a vector in two parts at (around) the middle part.
%
% Input argument:
% - vec: vector to be split up.
%
% Output
% - v1,v2: first and second half of vector
%
% Christian Utzerath 2017 (Donders Institute)


midindex =  floor (length (vec) / 2) ;
v1 = vec (1:midindex) ;
v2  = vec (midindex+1 :end) ;


end

