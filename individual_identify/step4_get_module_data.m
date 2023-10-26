%% extract the connectivity according to module
load('Community_detect.mat','com')

com_index = com.str_com;
com_num = length(unique(com.str_com));

roi_num = 146;
within_mask = ones(roi_num);
between_mask = diag(ones(roi_num,1));
all_mask = [within_mask,between_mask,between_mask,between_mask;...
    between_mask,within_mask,between_mask,between_mask;...
    between_mask,between_mask,within_mask,between_mask;...
    between_mask,between_mask,between_mask,within_mask];
all_mask(1:length(all_mask)+1:end) = 0;

sub_num = 114;
data_path = [pwd '\data'];


%% get data
for i = 1:com_num
    % within
    node_num = sum(com_index == i);
    mask = logical(tril(ones(node_num),-1));
    sub_data = zeros(sub_num,node_num*(node_num-1)/2);
    for isub = 1:sub_num
        clear SIM
        load([data_path '\JS_KSDENSITY_256_Signal_LR_sub_' num2str(isub,'%04d') '.mat'])
        SIM = SIM.*all_mask;
        data = SIM(com_index == i,com_index == i);
        sub_data(isub,:) = data(mask);
    end
    sub_data(:,mean(abs(sub_data)) == 0) = [];
    save(['Sub_matrix_data_module_' num2str(i) '_' num2str(i) '.mat'],'sub_data')
        
    % between
    for j = i+1:com_num
        node_num_i = sum(com_index == i);
        node_num_j = sum(com_index == j);
        sub_data = zeros(sub_num,node_num_j*node_num_i);
        for isub = 1:sub_num
            clear SIM
            load([data_path '\JS_KSDENSITY_256_Signal_LR_sub_' num2str(isub,'%04d') '.mat'])
            SIM = SIM.*all_mask;
            data = SIM(com_index == i,com_index == j);
            sub_data(isub,:) = data(:);
        end
        sub_data(:,mean(abs(sub_data)) == 0) = [];
        save(['Sub_matrix_data_module_' num2str(i) '_' num2str(j) '.mat'],'sub_data')
    end
end