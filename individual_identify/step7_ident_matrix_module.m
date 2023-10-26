%% identify the subjects
%% module connectivity
sub_num = 114;

load('Community_detect.mat','com')
com_num = length(unique(com.str_com));

for i = 1:com_num
    for j = i:com_num
        load(['Sub_matrix_data_module_' num2str(i) '_' num2str(i) '.mat'],'sub_data')
        [Acc_con, ~] = gretna_individual_identification(sub_data(1:sub_num/2,:),sub_data(1+sub_num/2:end,:),1,'Spearman');
        Results.Acc_con(i,j,1) = Acc_con.real;
        [Acc_con, ~] = gretna_individual_identification(sub_data(1+sub_num/2:end,:),sub_data(1:sub_num/2,:),1,'Spearman');
        Results.Acc_con(i,j,2) = Acc_con.real;
    end
end
Results.Acc_con = mean(Results.Acc_con,3);

save('Indi_matrix_module.mat','Results')