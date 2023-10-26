%% para
load('gene_a2009s_data.mat')
gene_data = corr(gene_data');
gene_data = (gene_data + gene_data')/2;
gene_data(1:length(gene_data)+1:end) = nan;
gene_data = gene_data(1:length(gene_data)/2,1:length(gene_data)/2);
gene_data = gene_data(logical(tril(ones(length(gene_data)),-1)));
del_edge = find(isnan(gene_data));
gene_data(del_edge) = [];

load('SC_data_a2009s.mat','Data_all')
modal = {'CT','FD','GI','SD'};
for i = 1:length(modal)
    Data.(modal{i}) = nanmean(Data_all.(modal{i}),3);
    Data.(modal{i}) = Data.(modal{i})(1:length(Data.(modal{i}))/2,1:length(Data.(modal{i}))/2);
end

load('spin_per_a2009s.mat')
perm_id = perm_id(1:size(perm_id,1)/2,:);

results = struct;

%% cal
for i = 1:length(modal)
    data = Data.(modal{i});
    data = data(logical(tril(ones(length(data)),-1)));
    data(del_edge) = [];
    
    % real
    results(i).modal = modal{i};
    results(i).data = data;
    results(i).gene_data = gene_data;
    [results(i).r,results(i).p] = corr(data,gene_data,'type','Spearman');
    
    % spin
    r_spin = zeros(size(perm_id,2),1);
    for i_per = 1:size(perm_id,2)
        if mod(i_per,floor(size(perm_id,2)/10))==0, disp(['permutation ' num2str(i_per) ' of ' num2str(size(perm_id,2)) ' ' modal{i} '  |' datestr(clock)]); end
        data_spin = Data.(modal{i});
        data_spin = data_spin(perm_id(:,i_per),perm_id(:,i_per));
        data_spin = data_spin(logical(tril(ones(length(data_spin)),-1)));
        data_spin(del_edge) = [];
        r_spin(i_per) = corr(data_spin,gene_data,'type','Spearman');
    end
    results(i).r_spin = r_spin;
    results(i).p_spin = (sum(abs(r_spin) >= abs(results(i).r))+1)/(size(perm_id,2)+1);
end

%% save
save('gene_all_pearson_co_edge_corr_left.mat','results')