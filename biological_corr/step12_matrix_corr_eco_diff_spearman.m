%% cal the feature contribution
load('EcoCorrMat43.mat')

load('corr_matrix_eco_spearman.mat')
results_ori = results;

modal = {'CT','FD','GI','SD'};

results = struct;

%% cal
for i = 1:length(modal)
    results(i).modal = modal{i};
    results(i).topo_name = topo_name;
    results(i).eco_r_ori = results_ori(i).eco_r;
    mask = results_ori(i).index;
    results(i).index = mask;
    
    % real
    data_JS = results_ori(i).data_JS;
    results(i).data_JS = data_JS;
    results(i).data_eco_ori = results_ori(i).data_eco;
    results(i).data_eco = zeros(length(data_JS),length(topo_name));
    
    for iname = 1:length(topo_name)
        data_eco = eco_data_z;
        data_eco(:,iname) = [];
        data_eco = corr(data_eco','row','pairwise','type','Spearman');
        data_eco = (data_eco + data_eco')/2;
        data_eco(1:length(data_eco)+1:end) = 0;
    
        results(i).data_eco(:,iname) = data_eco(mask);
    end
    
    results(i).eco_r = corr(data_JS,results(i).data_eco,'type','Spearman')';
    results(i).eco_r_diff = results(i).eco_r_ori - results(i).eco_r;
    
    [~,name_order] = sort(results(i).eco_r_diff,'descend');
    results(i).topo_name = [topo_name(name_order),num2cell(results(i).eco_r_diff(name_order)),...
        num2cell(results(i).eco_r_diff(name_order)/results(i).eco_r_ori)];
end

%% save
save('corr_matrix_eco_diff_spearman.mat','results')