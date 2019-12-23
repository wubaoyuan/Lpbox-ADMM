%% computes the GCO energy as defined by the GCO package
% ++ make sure x_sol starts with a label >= 1. 0 is not allowed as a label

function [energy] = compute_GCO_energy(unaryCosts,binary_cost_mat,x_sol)
[num_labels,num_pixels] = size(unaryCosts);
h = GCO_Create(num_pixels,num_labels);

% convert to int32
GCO_SetDataCost(h,int32(unaryCosts));

% set the smooth cost: Potts model
GCO_SetSmoothCost(h,ones(num_labels)-eye(num_labels));

% set the binary costs: cost of assigning neighboring pixels to different labels
GCO_SetNeighbors(h,binary_cost_mat);

% set the labeling
GCO_SetLabeling(h,x_sol);
[energy] = GCO_ComputeEnergy(h);

% Delete the GCoptimization object when finished
GCO_Delete(h);
return