function [condition_mat,sel_con_mat,condition_betas,condition_inds,baseline] = beta_matrix(nrun,ncon,nreg)
% Produces lists and matrixes of beta indices for your experiment.
% It is assumes that every session has the same conditions and noise
% regressors, and that all conditions are entered (i.e. do not skip
% conditions).
%
% Input:
% - nrun: number of sessions of experiment
% - ncon: number of conditions 
% - nreg: number of noise regressors (e.g., motion parameters)
%
% Output:
% - condition_mat:   ncon x nrun matrix with betas per run, with numbers
%                    reflecting position in SPM design
% - sel_con_mat:     ncon x nrun matrix, where each number denotes number of
%                    beta only in the selected subset, good for looping
%                    through runs of an experiment
% - condition_betas: binary vector highlighting selected betas
% - condition_inds:  just the numbers of the selected betas (a reshaped
%                    sel_con_mat), good for reading out experiment
% - baseline:        indices of baseline betas (whether or not you need
%                     them)
%
% Christian Utzerath 2015 (Donders Institute)


%% Do
condition_betas = repmat([ones(1,ncon) zeros(1,nreg)],1,nrun);              % Condition betas = 1, noise regressors = 0
condition_inds  = find(condition_betas == 1);                               % Indices of betas that reflect conditions
condition_mat   = reshape(condition_inds,ncon,nrun);                        % Reformat to conditions x runs matrix
baseline        = numel(condition_betas)+1:numel(condition_betas)+nrun;     % These would refer to baseline/offset/mean

selected = 1:numel(condition_inds);
sel_con_mat = reshape(selected,ncon,nrun);


end

