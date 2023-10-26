%% extract the connectivity according to module
load('Community_detect_single.mat','com')
modal = {'CT','FD','GI','SD'};

roi_num = 148;
within_mask = ones(roi_num);
all_mask = within_mask;
all_mask(1:length(all_mask)+1:end) = 0;

sub_num = 650;
path = pwd;
cd ..
data_path = [pwd '\data'];
cd(path)


%% get data
for imodal = 1:length(modal)
    com_index = com.(modal{imodal});
    com_num = length(unique(com.(modal{imodal})));
    for i = 1:com_num
        % within
        node_num = sum(com_index == i);
        mask = logical(tril(ones(node_num),-1));
        sub_data = zeros(sub_num,node_num*(node_num-1)/2);
        for isub = 1:sub_num
            clear SIM
            load([data_path '\JS_KSDENSITY_256_Signal_LR_sub_' num2str(isub,'%04d') '.mat'])
            SIM = SIM(1+roi_num*(imodal-1):roi_num*imodal,1+roi_num*(imodal-1):roi_num*imodal);
            SIM = SIM.*all_mask;
            data = SIM(com_index == i,com_index == i);
            sub_data(isub,:) = data(mask);
        end
        sub_data(:,mean(abs(sub_data)) == 0) = [];
        save(['Sub_matrix_data_module_' modal{imodal} '_' num2str(i) '_' num2str(i) '.mat'],'sub_data')
        
        % between
        for j = i+1:com_num
            node_num_i = sum(com_index == i);
            node_num_j = sum(com_index == j);
            sub_data = zeros(sub_num,node_num_j*node_num_i);
            for isub = 1:sub_num
                clear SIM
                load([data_path '\JS_KSDENSITY_256_Signal_LR_sub_' num2str(isub,'%04d') '.mat'])
                SIM = SIM(1+roi_num*(imodal-1):roi_num*imodal,1+roi_num*(imodal-1):roi_num*imodal);
                SIM = SIM.*all_mask;
                data = SIM(com_index == i,com_index == j);
                sub_data(isub,:) = data(:);
            end
            sub_data(:,mean(abs(sub_data)) == 0) = [];
            save(['Sub_matrix_data_module_' modal{imodal} '_' num2str(i) '_' num2str(j) '.mat'],'sub_data')
        end
    end
end