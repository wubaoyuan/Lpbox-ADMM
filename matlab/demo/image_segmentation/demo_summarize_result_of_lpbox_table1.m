clear
cd '/data1/wub/lpbox_admm/demo/image_segmentation/result'

pNorm_list = [0.5, 1, 2, 5, 10]; 
node_list = [1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6];
energy_mat = zeros(numel(node_list), numel(pNorm_list)); 
for n_id = 1:numel(node_list)
    n = node_list(n_id);
    for p_id = 1:numel(pNorm_list)
        pNorm = pNorm_list(p_id);
        file_list = dir(['*n_', num2str(n),'_pNorm_', num2str(pNorm), '*.mat']); 
        num_files = numel(file_list);
        energy_list = zeros(num_files, 1);
        for f = 1:numel(file_list)
           load(file_list(f).name); 
           energy_list(f) = result_struct_admm.energy;
        end
        energy_mat(n_id, p_id) = min(energy_list); 
    end
end

