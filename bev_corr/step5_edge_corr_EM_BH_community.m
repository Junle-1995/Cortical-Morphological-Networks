%% corr the connectivity at the community level and behav data via
%% EM model (https://github.com/RaphaelLiegeois/FC-Behavior/)
%% para
load('Community_detect.mat')
com_num = length(unique(com.str_com));

rand_time = 10000;
max_iter = 100;
tol = 1e-4;

path = pwd;
sub_cov = [];
load('HCP_unrelated_BH_data.mat')
BH_index = 1:6;

results = struct;


%% data
cd(path)
Predictors = cell(com_num);
for i = 1:com_num
    for j = i:com_num
        load(['Sub_matrix_data_module_' num2str(i) '_' num2str(j) '.mat'],'sub_data')
        Predictors{i,j} = zscore(sub_data);
    end
end


%% prediction
for iBH = BH_index
    results(iBH).BH_index = iBH;
    results(iBH).max_iter = max_iter;
    
    cd(path)
    bev_data = BH_unrelated_data(:,iBH);

    for i = 1:com_num
        for j = i:com_num
            [flag, m2, ~, ~, ~, ~, ~, ~] = get_individual_effects(bev_data, [], Predictors{i,j}, [], 0, tol, max_iter, 0);
            results(iBH).conver_flag(i,j) = flag;
            results(iBH).variance(i,j) = m2;
        end
    end
    
    for irand = 1:rand_time
        if mod(irand,floor(rand_time/10)) == 0
            disp(['Now calculating the data of randomization (' ...
                num2str(irand) '\' num2str(rand_time) ')  in ' num2str(iBH) ' |' datestr(clock)])
        end
        
        rand_index =randperm(length(bev_data));
        while isequal(rand_index,1:length(bev_data)), rand_index = randperm(length(bev_data)); end
        
        bev_data_rand = bev_data(rand_index);
        
        for i = 1:com_num
            for j = i:com_num
                [~, m2, ~, ~, ~, ~, ~, ~] = get_individual_effects(bev_data_rand, [], Predictors{i,j}, [], 0, tol, max_iter, 0);
                results(iBH).variance_rand(i,j,irand) = m2;
            end
        end
    end
    
    results(iBH).p = (sum(results(iBH).variance_rand >= repmat(results(iBH).variance,1,1,rand_time),3) + 1)/(rand_time + 1);
    
    cd(path)
    save('Edge_BH_EM_corr_community.mat','results')
end