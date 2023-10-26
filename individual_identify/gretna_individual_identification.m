function [Accuracy, P] = gretna_individual_identification(Data_target,Data_base,N_iter,Correlation)

%==========================================================================
% This function is used to identify individuals based on repeatedly observed
% data for the same set of participants.
%
%
% Syntax: [Accuracy, P] = gretna_individual_identification(Data_base,Data_target,N_iter)
%
% Inputs:
% Data_target:
%          N*P array with N being subjects and P being variables.
% Data_base:
%          N*P array with N being subjects and P being variables.
%    N_iter:
%          The number of iterations for nonparametric permutation test.
% Correlation(opt):
%          The correlative analysis methods.
%          'Pearson'£¨default£©
%          'Kendall'
%          'Spearman'
% Outputs:
% Accuracy.real:
%          The accuracy of real idetification.
% Accuracy.rand:
%          The accuracy of idetification for each iteraton.
%        P:
%          The significance of the real accuracy.
%
% Reference:
%          Emily S Finn et al.(2015): Functional connectome fingerprinting:
%          identifying individuals using patterns of brain connectivity.
%
% Zhen Li,     IBRR, SCNU, GuangZhou, 2020/09/23, zhen.li.dada@gmail.com
% Jinhui Wang, IBRR, SCNU, GuangZhou, 2019/12/26, jinhui.wang.1982@gmail.com
%==========================================================================

[N_sub_target,N_var_target] = size(Data_target);
[N_sub_base,N_var_base]     = size(Data_base);

if N_sub_base ~= N_sub_target
    error('The number of subjects is not equal between the two obervations!')
end

if N_var_base ~= N_var_target
    error('The number of variables is not equal between the two obervations!')
end

if nargin < 4
    R = corr(Data_target',Data_base','rows','pairwise');
else
    R = corr(Data_target',Data_base','rows','pairwise','type',Correlation);
end


Ind_ID_real = 1:N_sub_target;
[~,Ind_ID_predicted] = max(R);

Accuracy = struct;
Accuracy.real = sum(Ind_ID_real == Ind_ID_predicted)/N_sub_target;

% permutation
Accuracy.rand = zeros(N_iter,1);

for i_iter = 1:N_iter
    Ind_ID_rand = randperm(N_sub_target);
    Accuracy.rand(i_iter,1) = sum(Ind_ID_rand == Ind_ID_predicted)/N_sub_target;
end

P = (sum(Accuracy.rand > Accuracy.real))/(N_iter+1);

return