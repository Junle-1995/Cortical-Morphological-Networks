%% tidy the results of EM correlation
%% compare the results according to connectivity type and behavior type
load('Edge_community_BH_prediction_PCA_80_IQR.mat')

Results = struct;
for iBH = 1:6
    Results(iBH).BH_index = iBH;
end

module_num = max(cat(1,results(:).com_i));
for i = 1:length(results)
    for iBH = 1:6
        Results(iBH).R_mean(results(i).com_i,results(i).com_j) = ...
            results(i).R_mean(iBH);
        Results(iBH).p(results(i).com_i,results(i).com_j) = ...
            results(i).P_rand(iBH);
    end
end

for i = 1:length(Results)
    R_mean = [];
    R_mean(1:module_num,1) = diag(Results(i).R_mean);
    for imodule = 1:module_num
        R_mean(end+1:end+module_num-imodule) = Results(i).R_mean(imodule,imodule+1:end);
    end
    Results(i).R_mean = R_mean;
    
    p = [];
    p(1:module_num,1) = diag(Results(i).p);
    for imodule = 1:module_num
        p(end+1:end+module_num-imodule) = Results(i).p(imodule,imodule+1:end);
    end
    Results(i).p = p;
end

%% compare
Data = cat(1,Results(:).R_mean);
Inter = repmat([ones(module_num,1);2*ones(module_num*(module_num-1)/2,1)],length(Results),1);
Type = [];
for i = 1:length(Results)
    Type(end+1:(end+module_num+module_num*(module_num-1)/2),1) = i;
end
Data = [Data,Inter,Type];
out = SRH_test(Data,'Inter','Type');
P.Inter = out{1,5};
P.Type = out{2,5};
P.Interaction = out{3,5};

% type
P.Type_post = zeros(length(Results));
P.Type_ranksum = zeros(length(Results));
for itype = 1:length(Results)
    for jtype = itype+1:length(Results)
        [P.Type_post(itype,jtype),~,s] = ranksum(Results(itype).R_mean,Results(jtype).R_mean);
        P.Type_ranksum(itype,jtype) = s.zval;
        P.Type_ranksum(jtype,itype) = -s.zval;
    end
end
P.Type_post = P.Type_post + P.Type_post';
P.type_post_FDR = gretna_FDR(P.Type_post(logical(triu(ones(length(Results)),1))),0.05);

save('Edge_community_BH_prediction_PCA_80_IQR_compare.mat','P','Results')