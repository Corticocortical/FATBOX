function [iseven] = isEven(number)
%  This function returns a boolean indicating whether the input
% value is even (i.e. divisible by 2). 
%
% Output: 
% iseven: boolean
%
% Input: 
% number: scalar
%
% Christian Utzerath 2015 (Donders Institute)

if mod(number,2)
    iseven = 0;
else
    iseven = 1;
end

end

