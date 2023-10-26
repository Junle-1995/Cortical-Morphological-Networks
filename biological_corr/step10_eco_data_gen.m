%% gen eco data and matrix
load('scholtens2018neuroimage_economo.mat')
roi_num = 43;

topo_name = fieldnames(ECONOMO_43);

eco_data = zeros(roi_num,length(topo_name)-1);
for itopo = 1:length(topo_name)-1
    eco_data(:,itopo) = ECONOMO_43.(topo_name{itopo});
end

eco_data_z = (eco_data - repmat(nanmean(eco_data),roi_num,1))./(repmat(nanstd(eco_data),roi_num,1));

EcoCorrMat = corr(eco_data_z','type','Spearman','rows','pairwise');
EcoCorrMat = EcoCorrMat + EcoCorrMat';
EcoCorrMat = EcoCorrMat/2;
EcoCorrMat(1:length(EcoCorrMat)+1:end) = 0;
topo_name(end) = [];

save('EcoCorrMat43.mat','EcoCorrMat','eco_data_z','eco_data','topo_name')