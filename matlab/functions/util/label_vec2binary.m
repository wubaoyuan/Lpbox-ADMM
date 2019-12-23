%% converts the label_vec into a binary matrix

function [bin_label_array] = label_vec2binary(label_vec,num_classes)
num_nodes = numel(label_vec);
bin_label_array = zeros(num_nodes,num_classes);

for c = 1:num_classes
    bin_label_array(label_vec==c,c) = 1;
end
return
