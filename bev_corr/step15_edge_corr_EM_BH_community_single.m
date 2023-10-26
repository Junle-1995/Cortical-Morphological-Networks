%% corr the connectivity at the community level and behav data via
%% EM model (https://github.com/RaphaelLiegeois/FC-Behavior/)
%% para
load('Community_detect_single.mat')
modal = {'CT','FD','GI','SD'};

max_iter = 100;
tol = 1e-4;

path = pwd;
sub_cov = [];
load('HCP_unrelated_BH_data.mat')
BH_index = 1:6;

results = struct;


%% cal
for imodal = 1:length(modal)
    
    % data
    cd(path)
    com_num = length(unique(com.(modal{imodal})));
    Predictors = cell(com_num);
    for i = 1:com_num
        for j = i:com_num
            load(['Sub_matrix_data_module_' modal{imodal} '_' num2str(i) '_' num2str(j) '.mat'],'sub_data')
            Predictors{i,j} = zscore(sub_data);
        end
    end
    
    
    % prediction
    for iBH = BH_index
        results.(modal{imodal})(iBH).BH_index = iBH;
        results.(modal{imodal})(iBH).max_iter = max_iter;
        
        cd(path)
        bev_data = BH_unrelated_data(:,iBH);
        
        for i = 1:com_num
            for j = i:com_num
                [flag, m2, ~, ~, ~, ~, ~, ~] = get_individual_effects(bev_data, [], Predictors{i,j}, [], 0, tol, max_iter, 0);
                results.(modal{imodal})(iBH).conver_flag(i,j) = flag;
                results.(modal{imodal})(iBH).variance(i,j) = m2;
            end
        end
    end
end

cd(path)
save('Edge_BH_EM_corr_community_single.mat','results')