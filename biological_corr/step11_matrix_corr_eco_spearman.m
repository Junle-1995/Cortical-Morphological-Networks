%% cal correlation between JS matrix and eco matrix
JS_data = load('con_unrelated_mean_JS.mat');
load('EcoCorrMat43.mat')

modal = {'CT','FD','GI','SD'};
load('spin_per_a2009s.mat')

load('ECO_Left_a2009s.mat')
EC_label = LABEL_MAPPING;
load('ECO_Right_a2009s.mat')
EC_label = [EC_label; LABEL_MAPPING];
EC_label = EC_label.economo_id;
for i = 1:length(EC_label)
    if isempty(EC_label{i}), EC_label{i} = 0; end
end
EC_label = cat(1,EC_label{:});
mask = logical(triu(ones(43),1));

per_time = 10000;

results = struct;

%% cal
for i = 1:length(modal)
    results(i).modal = modal{i};
    
    % real
    data_eco = EcoCorrMat;
    data_JS_ori = JS_data.([modal{i} '_mean_JS']);
    data_JS = zeros(length(data_eco));
    for ieco = 1:length(data_JS)
        data_i = data_JS_ori(EC_label == ieco,EC_label == ieco);
        data_i(1:length(data_i)+1:end) = nan;
        if length(data_i) ~= 1
            data_JS(ieco,ieco) = nanmean(data_i(:));
        end
        
        for jeco = ieco+1:length(data_JS)
            data_i = data_JS_ori(EC_label == ieco,EC_label == jeco);
            if ~isempty(data_i)
                data_JS(ieco,jeco) = nanmean(data_i(:));
                data_JS(jeco,ieco) = nanmean(data_i(:));
            end
        end
    end
    
    results(i).data_JS = data_JS(mask & data_JS~= 0);
    results(i).data_eco = data_eco(mask & data_JS~= 0);
    results(i).index = mask & data_JS~= 0;
    
    results(i).eco_r = corr(data_JS(mask & data_JS~= 0),data_eco(mask & data_JS~= 0),'type','Spearman');
    
    % spin
    results(i).eco_r_spin = zeros(per_time,1);
    for i_per = 1:per_time
        if mod(i_per,per_time/10) == 0, disp(['Now calculating the data in ' modal{i} ' (' num2str(i_per) '\' num2str(per_time) ')  |' datestr(clock)]); end
        
        data_JS_ori_spin = data_JS_ori(perm_id(:,i_per),perm_id(:,i_per));
        data_JS_spin = zeros(length(data_eco));
        for ieco = 1:length(data_JS_spin)
            data_i = data_JS_ori_spin(EC_label == ieco,EC_label == ieco);
            data_i(1:length(data_i)+1:end) = nan;
            if length(data_i) ~= 1
                data_JS_spin(ieco,ieco) = nanmean(data_i(:));
            end
            
            for jeco = ieco+1:length(data_JS_spin)
                data_i = data_JS_ori_spin(EC_label == ieco,EC_label == jeco);
                if ~isempty(data_i)
                    data_JS_spin(ieco,jeco) = nanmean(data_i(:));
                    data_JS_spin(jeco,ieco) = nanmean(data_i(:));
                end
            end
        end
        
        results(i).eco_r_spin(i_per,1) = corr(data_JS_spin(mask & data_JS_spin ~= 0),data_eco(mask & data_JS_spin ~= 0),'type','Spearman');
    end
    results(i).eco_p = (sum(abs(results(i).eco_r_spin) >= abs(results(i).eco_r))+1)/(per_time+1);
end

%% save
save('corr_matrix_eco_spearman.mat','results')