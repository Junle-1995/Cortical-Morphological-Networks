%% compare the results of within-module and between-module
load('Indi_matrix_module.mat')

between_mask = logical(triu(ones(length(Results.Acc_con)),1));

[p,h,stats] = ranksum(diag(Results.Acc_con),Results.Acc_con(between_mask),'method','approximate');

Results.p = p;
Results.z = stats.zval;

save('Indi_matrix_module.mat','Results')