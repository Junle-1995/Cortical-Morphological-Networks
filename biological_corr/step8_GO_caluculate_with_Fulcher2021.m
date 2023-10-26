%% para
%% required toolbox: https://github.com/benfulcher/GeneCategoryEnrichmentAnalysis
load('gene_all_pearson_co_edge_corr_left_gene_null.mat','results')
Index = 1;

sizeFilter = [10 1000];
load('gene_name_ID.mat')
geneEntrezIDs = gene_name_ID(~isnan(gene_name_ID));


%% get GO table
GOTable = GetFilteredGOData('human-direct','biological_process',sizeFilter,geneEntrezIDs);


%% calculate the real GO score
GO_results = struct;
geneScores = results(Index).r_change;
geneScores = geneScores(~isnan(gene_name_ID));

disp(['Now calculating the real GO score   |' datestr(clock)])
for iGO = 1:size(GOTable,1)
    GO_results(iGO).GOlabel = GOTable.GOIDlabel{iGO};
    GO_results(iGO).GOname = GOTable.GOName{iGO};
    GO_results(iGO).GOsize = GOTable.size(iGO);
    
    GO_results(iGO).GOscore = mean(geneScores(ismember(geneEntrezIDs,GOTable.annotations{iGO})));
end


%% get the permutated GO score
per_time = size(results(Index).gene_per_r_change,2);
for iper = 1:per_time
    disp(['Now calculating the random GO score (' num2str(iper) '/' num2str(per_time) ')  |' datestr(clock)])
    
    geneScores = results(Index).gene_per_r_change(:,iper);
    geneScores = geneScores(~isnan(gene_name_ID));
    
    for iGO = 1:size(GOTable,1)
        GO_results(iGO).GOscore_rand(iper,1) = mean(geneScores(ismember(geneEntrezIDs,GOTable.annotations{iGO})));
    end
end


%% calculated the p-value of each GO
for iGO = 1:size(GOTable,1)
    GO_results(iGO).GOscore_per_p = mean(GO_results(iGO).GOscore_rand >= GO_results(iGO).GOscore);
end

per_p_fdr = mafdr(cat(2,GO_results(:).GOscore_per_p),'BHFDR',true,'showPlot',false);
[~,change_index_per] = sort(per_p_fdr,'ascend');

for iGO = 1:size(GOTable,1)
    GO_results(iGO).GOscore_per_p_FDR = per_p_fdr(iGO);
end
GO_results_per = GO_results(change_index_per);


%% save
save('gene_all_pearson_co_edge_corr_left_gene_null_GOterm.mat','GO_results_per')