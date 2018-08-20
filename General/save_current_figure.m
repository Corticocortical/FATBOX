function [ output_args ] = save_current_figure(folder,name,trans)
% Shorthand to save a figure using hg_export.
% hg_export must be on path - might not be by default.
% Uses PNG.
%
% Input arguments:
%  - folder: string. Where to save figure.
%  - name: string. How to name file. No file format.
%  - trans: 0/1. Transparent or not.
%
% Christian Utzerath 2015 (Donders Institute)


%% Input args
returndir = pwd;
cd(folder)

if trans
    t = '-transparent';
else
    t = '';
end


%% Assemble command
%cellfun(@export_fig,{name t}) %  weird
command = ['export_fig ''' name ''' ' t];
eval(command)

cd(returndir)


end

