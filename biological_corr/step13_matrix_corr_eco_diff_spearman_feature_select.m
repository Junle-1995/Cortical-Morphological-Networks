%% select the feature significantly contribute to the correlation
file_name = 'corr_matrix_eco_diff_spearman.mat';

load(file_name)

% CT
diff_z = zscore(cat(1,results(1).topo_name{:,2}));
results(1).topo_name(:,4) = num2cell(diff_z);
eco_r_diff_p = 1 - cdf('Normal',cat(1,results(1).topo_name{:,4}),0,1);
results(1).topo_name(:,5) = num2cell(eco_r_diff_p);

% save
save(file_name,'results')