%% cal the feature contribution
JS_data = load('con_unrelated_mean_JS.mat');
load('neuro_trans_net_spearman.mat')
TM_data = TM_data.a2009s;

load('corr_matrix_spearman.mat')
results_ori = results;

modal = {'CT','FD','GI','SD'};
load('spin_per_a2009s.mat')
mask = logical(triu(ones(148),1));

per_time = 10000;

results = struct;

%% cal
for i = 1:length(modal)
    results(i).modal = modal{i};
    results(i).trans_name = trans_name;
    results(i).r_ori = results_ori(i).trans_r;
    
    % real
    data_JS = JS_data.([modal{i} '_mean_JS']);
    results(i).data_JS = data_JS(mask);
    results(i).data_trans = zeros(length(data_JS(mask)),length(trans_name));
    
    for iname = 1:length(trans_name)
        data_trans = TM_data;
        data_trans(:,iname) = [];
        data_trans = ...
            (data_trans-repmat(nanmean(data_trans),size(data_trans,1),1))./repmat(nanstd(data_trans),size(data_trans,1),1);
        data_trans = corr(data_trans','row','pairwise','type','Spearman');
        data_trans = (data_trans + data_trans')/2;
        data_trans(1:length(data_trans)+1:end) = 0;
    
        results(i).data_trans(:,iname) = data_trans(mask);
    end
    
    results(i).trans_r = corr(data_JS(mask),results(i).data_trans,'type','Spearman')';
    results(i).trans_r_diff = results(i).r_ori - results(i).trans_r;
    
    [~,name_order] = sort(results(i).trans_r_diff,'descend');
    results(i).trans_name = [trans_name(name_order),num2cell(results(i).trans_r_diff(name_order)),...
        num2cell(results(i).trans_r_diff(name_order)/results(i).r_ori)];
end

%% save
save('corr_matrix_trans_diff_spearman.mat','results')