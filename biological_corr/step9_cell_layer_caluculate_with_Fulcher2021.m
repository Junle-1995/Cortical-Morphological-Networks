%% para
load('gene_all_pearson_co_edge_corr_left_gene_null.mat','results')
Index = 1;

per_time = size(results(Index).gene_per_r_change,2);


%% get type
load('gene_cell_type.mat')
Type_assign = cell_type_assign;
load('gene_layer.mat')
Type_assign.layer13 = gene_assign.layer13;
Type_assign.layer4 = gene_assign.layer4;
Type_assign.layer56 = gene_assign.layer56;

Type_name = fieldnames(Type_assign);


%% calculate the real score
GO_results = struct;
geneScores = results(Index).r_change;

disp(['Now calculating the real score   |' datestr(clock)])
for iGO = 1:length(Type_name)
    GO_results(iGO).GOname = Type_name{iGO};
    GO_results(iGO).GOsize = length(Type_assign.(Type_name{iGO}));
    
    GO_results(iGO).GOscore = mean(geneScores(Type_assign.(Type_name{iGO})));
end


%% get the permutated score
for iper = 1:per_time
    disp(['Now calculating the random score (' num2str(iper) '/' num2str(per_time) ')  |' datestr(clock)])
    
    geneScores = results(Index).gene_per_r_change(:,iper);
    
    for iGO = 1:length(Type_name)
        GO_results(iGO).GOscore_rand(iper,1) = mean(geneScores(Type_assign.(Type_name{iGO})));
    end
end


%% calculated the p-value of each GO
for iGO = 1:length(GO_results)
    GO_results(iGO).GOscore_per_p = mean(GO_results(iGO).GOscore_rand >= GO_results(iGO).GOscore);
end


%% save
save('gene_all_pearson_co_edge_corr_left_gene_null_cell_layer.mat','GO_results')