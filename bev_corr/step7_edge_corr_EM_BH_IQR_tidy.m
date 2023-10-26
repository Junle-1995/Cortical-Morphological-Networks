%% tidy the results of EM correlation
%% check the conver_flag
%% compare the results according to connectivity type and behavior type
load('Edge_BH_EM_corr_community_IQR.mat')

for i = 1:length(results)
    module_num = length(results(i).conver_flag);
    
    con_flag = [];
    con_flag(1:module_num,1) = diag(results(i).conver_flag);
    for imodule = 1:module_num
        con_flag(end+1:end+module_num-imodule) = results(i).conver_flag(imodule,imodule+1:end);
    end
    results(i).conver_flag = con_flag;
    if sum(results(i).conver_flag == 0) > 1, error('No convert exist!'), end
    
    var = [];
    var(1:module_num,1) = diag(results(i).variance);
    for imodule = 1:module_num
        var(end+1:end+module_num-imodule) = results(i).variance(imodule,imodule+1:end);
    end
    results(i).variance = var;
    
    p = [];
    p(1:module_num,1) = diag(results(i).p);
    for imodule = 1:module_num
        p(end+1:end+module_num-imodule) = results(i).p(imodule,imodule+1:end);
    end
    results(i).p = p;
end

%% compare
Data = cat(1,results(:).variance);
Inter = repmat([ones(module_num,1);2*ones(module_num*(module_num-1)/2,1)],length(results),1);
Type = [];
for i = 1:length(results)
    Type(end+1:(end+module_num+module_num*(module_num-1)/2),1) = i;
end
Data = [Data,Inter,Type];
out = SRH_test(Data,'Inter','Type');
P.Inter = out{1,5};
P.Type = out{2,5};
P.Interaction = out{3,5};

% type
P.Type_post = zeros(length(results));
P.Type_ranksum = zeros(length(results));
for itype = 1:length(results)
    for jtype = itype+1:length(results)
        [P.Type_post(itype,jtype),~,s] = ranksum(results(itype).variance,results(jtype).variance);
        P.Type_ranksum(itype,jtype) = s.zval;
        P.Type_ranksum(jtype,itype) = -s.zval;
    end
end
P.Type_post = P.Type_post + P.Type_post';
P.type_post_FDR = gretna_FDR(P.Type_post(logical(triu(ones(length(results)),1))),0.05);

save('Edge_BH_EM_corr_community_compare_IQR.mat','P','results')