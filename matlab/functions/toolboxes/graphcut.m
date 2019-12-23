function [E_a, labels_a, E_swap, labels_swap] = graphcut(unaryCosts, weights)
%%%%%%%%%%%%%%%%%%%%
% MAP estimate of the CRF cost using the graph-cut algorithm from:
% Multi-label optimization http://vision.csd.uwo.ca/code/ 
% an implementation of alpha-expansion and alpha-beta
% 
% Angjoo Kanazawa 5/3/'12
%%%%%%%%%%%%%%%%%%%%
    L = 2; % binary label for foreground background segmentation
    N = size(weights, 1);
    h = GCO_Create(N, L); 
    % set GCO parameters:
    GCO_SetDataCost(h, unaryCosts);
    % GCO_SetSmoothCost use default Potts model
    GCO_SetNeighbors(h, weights);
    GCO_Expansion(h);
    labels_a = GCO_GetLabeling(h);
    [E_a D S l] = GCO_ComputeEnergy(h);
    GCO_Delete(h);
    
    h2 = GCO_Create(N, L); 
    % set GCO parameters:
    GCO_SetDataCost(h2, unaryCosts);
    % GCO_SetSmoothCost use default Potts model
    GCO_SetNeighbors(h2, weights);
    GCO_Swap(h2);
    labels_swap = GCO_GetLabeling(h2);
    [E_swap D S l] = GCO_ComputeEnergy(h2);
    GCO_Delete(h2);
    labels_a=double(labels_a)-1;
    labels_swap=double(labels_swap)-1;
end
