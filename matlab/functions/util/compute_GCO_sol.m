function [GCO_label_vec,E_GCO] = compute_GCO_sol(unaryCosts,binary_cost_mat)
[num_labels,num_pixels] = size(unaryCosts);
h = GCO_Create(num_pixels,num_labels);

% convert to int32
GCO_SetDataCost(h,int32(unaryCosts));

% set the smooth cost: Potts model
GCO_SetSmoothCost(h,ones(num_labels)-eye(num_labels));

% set the binary costs: cost of assigning neighboring pixels to different labels
GCO_SetNeighbors(h,binary_cost_mat);

% compute optimal labeling via alpha-beta expansion 
GCO_Expansion(h);

% GCO_Swap(h);
GCO_label_vec = double(GCO_GetLabeling(h));
[E_GCO] = GCO_ComputeEnergy(h);

% Delete the GCoptimization object when finished
GCO_Delete(h);
return;
