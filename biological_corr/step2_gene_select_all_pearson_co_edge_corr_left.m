%% para
load('gene_all_pearson_co_edge_corr_left.mat')

load('gene_a2009s_data.mat')

load('SC_data_a2009s.mat','Data_all')
modal = {'CT','FD','GI','SD'};
for i = 1:length(modal)
    Data.(modal{i}) = nanmean(Data_all.(modal{i}),3);
    Data.(modal{i}) = Data.(modal{i})(1:length(Data.(modal{i}))/2,1:length(Data.(modal{i}))/2);
end
modal_index = 1;

z_thr = icdf('normal',0.95,0,1);

%% cal
for i = modal_index
    data = results(i).data;
    
    % gene
    r_gene_i = zeros(length(gene_name),1);
    for i_gene = 1:length(gene_name)
        if mod(i_gene,floor(length(gene_name)/10))==0, disp(['cal ' num2str(i_gene) ' of ' num2str(length(gene_name)) ' ' modal{i} '  |' datestr(clock)]); end
        
        gene_data_i = gene_data;
        gene_data_i(:,i_gene) = [];
        gene_data_i = corr(gene_data_i');
        gene_data_i = (gene_data_i + gene_data_i')/2;
        gene_data_i(1:length(gene_data_i)+1:end) = nan;
        gene_data_i = gene_data_i(1:length(gene_data_i)/2,1:length(gene_data_i)/2);
        gene_data_i = gene_data_i(logical(tril(ones(length(gene_data_i)),-1)));
        del_edge = find(isnan(gene_data_i));
        gene_data_i(del_edge) = [];
        
        r_gene_i(i_gene) = corr(data,gene_data_i,'type','spearman');
    end
    results(i).r_gene_i = r_gene_i;
    results(i).r_change = results(i).r - r_gene_i;
    
    [~,gene_index] = sort(results(i).r_change,'descend');
    results(i).gene_sort = gene_name(gene_index);
    
    gene_select = zscore(results(i).r_change);
    results(i).gene_select = gene_name(gene_select >= z_thr);
end

%% save
save('gene_all_pearson_co_edge_corr_left_gene.mat','results')