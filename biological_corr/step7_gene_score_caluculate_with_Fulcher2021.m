%% para
%% calculate the gene contribution in each permutation preserving spatial autocorrelation
load('gene_all_pearson_co_edge_corr_left_gene.mat')

load('gene_a2009s_data.mat')

load('SC_data_a2009s.mat','Data_all')
modal = {'CT','FD','GI','SD'};
for i = 1:length(modal)
    Data.(modal{i}) = nanmean(Data_all.(modal{i}),3);
    Data.(modal{i}) = Data.(modal{i})(1:length(Data.(modal{i}))/2,1:length(Data.(modal{i}))/2);
end
modal_index = 1;

load('spin_per_a2009s.mat')
perm_id(size(perm_id,1)/2+1:end,:) = [];
per_time = 10000;

%% cal
for i = modal_index
    data = Data.(results(i).modal);
    
    % all gene permutated r
    results(i).gene_all_per_r = zeros(per_time,1);
    gene_data_all = gene_data;
    gene_data_all = corr(gene_data_all');
    gene_data_all = (gene_data_all + gene_data_all')/2;
    gene_data_all(1:length(gene_data_all)+1:end) = nan;
    gene_data_all = gene_data_all(1:length(gene_data_all)/2,1:length(gene_data_all)/2);
    gene_data_all = gene_data_all(logical(tril(ones(length(gene_data_all)),-1)));
    
    for i_per = 1:per_time
        data_per = data(perm_id(:,i_per),perm_id(:,i_per));
        data_per = data_per(logical(tril(ones(length(data_per)),-1)));
        results(i).gene_all_per_r(i_per,1) = corr(gene_data_all,data_per,'type','Spearman','rows','pairwise');
    end
    
    % gene
    for i_gene = 1:length(gene_name)
        disp(['cal ' num2str(i_gene) ' of ' num2str(length(gene_name)) ' ' modal{i} '  |' datestr(clock)])
        
        r_gene_i = zeros(per_time,1);
        
        gene_data_i = gene_data;
        gene_data_i(:,i_gene) = [];
        gene_data_i = corr(gene_data_i');
        gene_data_i = (gene_data_i + gene_data_i')/2;
        gene_data_i(1:length(gene_data_i)+1:end) = nan;
        gene_data_i = gene_data_i(1:length(gene_data_i)/2,1:length(gene_data_i)/2);
        gene_data_i = gene_data_i(logical(tril(ones(length(gene_data_i)),-1)));
        
        for i_per = 1:per_time
            data_per = data(perm_id(:,i_per),perm_id(:,i_per));
            data_per = data_per(logical(tril(ones(length(data_per)),-1)));
            r_gene_i(i_per,1) = corr(data_per,gene_data_i,'type','Spearman','rows','pairwise');
        end
        results(i).gene_per_r_change(i_gene,:) = (results(i).gene_all_per_r - r_gene_i)';
    end
end

%% save
save('gene_all_pearson_co_edge_corr_left_gene_null.mat','results')