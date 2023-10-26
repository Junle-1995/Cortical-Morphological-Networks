%% corr the connectivity and behav data via gretna_prediction
%% para
load('Community_detect.mat')
com_num = length(unique(com.str_com));

K_fold = 10;
C_type = 'PCA';
C_thr = 80;
R_type = 'glm';
Ncomp = [];
Lambda = [];

rep_time = 100;
rand_time = 10000;

path = pwd;
load('HCP_unrelated_BH_data.mat')

results = struct;
results(1) = [];


%% prediction
Responses = BH_unrelated_data;
for i = 1:com_num
    for j = i:com_num
        results(end+1).com_i = i;
        results(end).com_j = j;
        
        cd(path)
        load(['Sub_matrix_data_module_' num2str(i) '_' num2str(j) '.mat'],'sub_data')
        Predictors = zscore(sub_data);
        
        for i_rep = 1:rep_time
            disp(['Now calculating the data of repeation (' ...
                num2str(i_rep) '\' num2str(rep_time) ')  in ' num2str(i) '_' num2str(j) ' |' datestr(clock)])
            
            Results = gretna_prediction_all_var_for_pca(Predictors, Responses, K_fold, C_type, C_thr, R_type, Ncomp, Lambda);
            results(end).Explained_variance(i_rep,1) = mean(Results.Explained_variance);
            results(end).R(:,i_rep) = Results.R;
            results(end).P(:,i_rep) = Results.P;
            results(end).Response_predicted(:,:,i_rep) = Results.Response_predicted;
        end
        results(end).MAE = mean(mean(abs(results(end).Response_predicted - repmat(Responses,1,1,rep_time)),3),1)';
        
        cd(path)
        results(end).MAE_rand = zeros(size(Responses,2),rand_time);
        for irand = 1:rand_time
            %         if mod(irand,floor(rand_time/10)) == 0
            disp(['Now calculating the data of randomization (' ...
                num2str(irand) '\' num2str(rand_time) ')  in ' num2str(i) '_' num2str(j) ' |' datestr(clock)])
            %         end
            
            rand_index = randperm(size(Responses,1));
            while isequal(rand_index,1:size(Responses,1)), rand_index =randperm(size(Responses,1)); end
            
            Responses_rand = Responses(rand_index,:);
            
            Results = gretna_prediction_all_var_for_pca(Predictors, Responses_rand, K_fold, C_type, C_thr, R_type, Ncomp, Lambda);
            results(end).R_rand(:,irand) = Results.R;
            results(end).MAE_rand(:,irand) = mean(abs(Results.Response_predicted - Responses_rand))';
            
            if mod(irand,floor(rand_time/10)) == 0, save(['Edge_community_BH_prediction_PCA_' num2str(C_thr) '_MAE.mat'],'results'), end
        end
        
        
        results(end).R_mean = mean(results(end).R,2);
        results(end).P_rand = (sum(results(end).R_rand >= repmat(results(end).R_mean,1,rand_time),2) + 1)/(rand_time + 1);
        results(end).P_MAE = (sum(results(end).MAE_rand <= repmat(results(end).MAE,1,rand_time),2) + 1)/(rand_time + 1);
        
        save(['Edge_community_BH_prediction_PCA_' num2str(C_thr) '_MAE.mat'],'results')
    end
end