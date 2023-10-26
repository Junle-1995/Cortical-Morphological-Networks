%% cal correlation between JS matrix and neurontrans matrix
JS_data = load('con_unrelated_mean_JS.mat');
load('neuro_trans_net_spearman.mat')

modal = {'CT','FD','GI','SD'};
load('spin_per_a2009s.mat')
mask = logical(triu(ones(148),1));

per_time = 10000;

results = struct;

%% cal
for i = 1:length(modal)
    results(i).modal = modal{i};
    
    % real
    data_JS = JS_data.([modal{i} '_mean_JS']);
    data_trans = TM_net.a2009s;
    
    results(i).data_JS = data_JS(mask);
    results(i).data_trans = data_trans(mask);
    
    results(i).trans_r = corr(data_JS(mask),data_trans(mask),'type','Spearman');
    
    % spin
    results(i).trans_r_spin = zeros(per_time,1);
    for i_per = 1:per_time
        data_JS_spin = data_JS(perm_id(:,i_per),perm_id(:,i_per));
        
        results(i).trans_r_spin(i_per,1) = corr(data_JS_spin(mask),data_trans(mask),'type','Spearman');
    end
    results(i).trans_p = (sum(abs(results(i).trans_r_spin) >= abs(results(i).trans_r))+1)/(per_time+1);
end

%% save
save('corr_matrix_spearman.mat','results')