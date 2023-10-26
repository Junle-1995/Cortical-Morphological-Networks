%% corr the connectivity and behav data via gretna_prediction
%% para
load('Community_detect_single.mat')
modal = {'CT','FD','GI','SD'};

K_fold = 10;
C_type = 'PCA';
C_thr = 80;
R_type = 'glm';
Ncomp = [];
Lambda = [];

rep_time = 100;

path = pwd;
load('HCP_unrelated_BH_data.mat')

results = struct;


%% prediction
for imodal = 1:length(modal)
    com_num = length(unique(com.(modal{imodal})));
    Responses = BH_unrelated_data;
    results.(modal{imodal}) = struct;
    results.(modal{imodal})(1) = [];
    
    for i = 1:com_num
        for j = i:com_num
            disp(['Now calculating the data in ' modal{imodal} '_' num2str(i) '_' num2str(j) ' |' datestr(clock)])
            
            results.(modal{imodal})(end+1).com_i = i;
            results.(modal{imodal})(end).com_j = j;
            
            cd(path)
            load(['Sub_matrix_data_module_' modal{imodal} '_' num2str(i) '_' num2str(j) '.mat'],'sub_data')
            Predictors = zscore(sub_data);
            
            for i_rep = 1:rep_time
                Results = gretna_prediction_all_var_for_pca(Predictors, Responses, K_fold, C_type, C_thr, R_type, Ncomp, Lambda);
                results.(modal{imodal})(end).Explained_variance(i_rep,1) = mean(Results.Explained_variance);
                results.(modal{imodal})(end).R(:,i_rep) = Results.R;
            end
            results.(modal{imodal})(end).R_mean = mean(results.(modal{imodal})(end).R,2);
        end
    end
end

cd(path)
save(['Edge_community_BH_prediction_PCA_' num2str(C_thr) '_single.mat'],'results')