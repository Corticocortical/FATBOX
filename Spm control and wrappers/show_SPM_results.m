function [ output_args ] = show_SPM_results(varargin)
% A shorthand to open the SPM results screen. Saves
% a lot of clicking if you want to do a lot of things
% through the GUI.
%
% Input argument(s):
% Required: 
% - folder. string. folder in which SPM.mat is stored. 
%   Optional:
%    - correction. the p correction to apply. If not supplied, will
%      default to p < 0.001 uncorrected with k = 0. Fields: 
%      * method: 'fwe' or 'none'
%      * p: p-value for correction
%      * k: voxel extent
%    - contrast number: the number of the contrast to open.
%
% Example usage: 
% show_SPM_results('myfolder/GLM')
% OR
% correction.method = 'none'; correction.p = 0.001; correction.k = 30;
% show_SPM_results('myfolder',correction)
%
% Christian Utzerath 2017 (Donders Institute)


%% Parse input
rtdir = pwd;
folder = varargin{1};
spmpath = [folder,filesep,'SPM.mat'];

% What correction to apply?
if numel(varargin) > 1
    correction = varargin{2};
else
    correction.method = 'none';
    correction.p      = 0.001;
    correction.k      = 0;
end

% What contrast to open first?
if numel(varargin) > 2
    cnoi = varargin{3}
else
    cnoi = 1;
end

%% Call SPM
clear jobs
jobs{1}.stats{1}.results.spmmat                     = cellstr(spmpath);
jobs{1}.stats{1}.results.conspec(1).titlestr        = ['Results view'];
jobs{1}.stats{1}.results.conspec(1).contrasts       = cnoi;              % contrast number
jobs{1}.stats{1}.results.conspec(1).threshdesc      = correction.method; % correction level
jobs{1}.stats{1}.results.conspec(1).thresh          = correction.p;      % p value threshold
jobs{1}.stats{1}.results.conspec(1).extent          = correction.k;      % extent threshold
jobs{1}.stats{1}.results.print                      = 0;                 % don't print (this should speed things up a bit)

spm('defaults','FMRI');
spm_jobman('initcfg');
spm_jobman('run',jobs);
cd(rtdir)

end

