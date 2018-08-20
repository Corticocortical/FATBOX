function [ram_consumed] = ram_consumption(ram)
% Adds up all the variables in the workspace to give an estimate of
% consumed RAM. Note that this function only reports what is in the
% accessible global workspace.
%
%
% Input:
% - ram. Output of the matlab whos function. Program cannot by itself call
%        this function.
%
% Output: 
% - ram_consumed: consumed RAM (in GB)
%
% Christian Utzerath 2016 (Donders Institute)
%%

weight = [];
for r = 1:numel(ram)
    weight(r) = ram(r).bytes;
end

ram_consumed = sum(weight)/1024/1024/1024;

end

