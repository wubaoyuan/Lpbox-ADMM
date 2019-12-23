clear
cd '/data1/wub/lpbox_admm/demo/clustering/results'
dataset_name = 'glass';
cd dataset_name
pNorm_list = [0.5, 1, 2, 5, 10];
record_matrix = zeros(numel(pNorm_list),4);
for p_id = 1:numel(pNorm_list)
	pNorm = pNorm_list(p_id);
	file_list = dir(['taskID_*', dataset_name, '*_Lp_', num2str(pNorm), '*.mat']);
	num_files = numel(file_list);
	obj_list = zeros(num_files, 1);
        ri_list = zeros(num_files, 1);
        record_list = zeros(num_files, 4);
	for f = 1:numel(file_list)
	   load(file_list(f).name);
           record_list(f,:) = resultStruct.record_list;
	   obj_list(f) = resultStruct.record_list(3);
	   ri_list(f) = resultStruct.record_list(1);
        end
	[~, max_index] = min(obj_list);
        record_matrix(p_id,:) = record_list(max_index,:);
end
