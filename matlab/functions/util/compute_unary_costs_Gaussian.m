%% computes unary costs based on how different the colors are to a Gaussian distribution
% gray scale images for now

function [gCosts] = compute_unary_costs_Gaussian(I,gMeans,gStd)
num_pixels = numel(I);
num_classes = numel(gStd);

gCosts = zeros(num_classes,num_pixels);
for c = 1:num_classes
    gCosts(c,:) = 0.5*log(2*pi)+log(gStd(c))+0.5*((I(:)-gMeans(c)).^2)/(gStd(c)^2);
end

return
