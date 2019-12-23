function [binary_cost_mat,time_elapsed] = generate_binary_cost(I,params)
if ~exist('params','var') || isempty(params)
    params = struct('sigma',std(I(:)),'neighborhood_edge',1);
end

tic; 

[nrows,ncols,num_colors] = size(I);
[pairs] = generate_pixel_pairs(nrows,ncols,params.neighborhood_edge);

% compute the feature differences
if num_colors==1
    diff_vec = ((I(pairs(:,1))-I(pairs(:,2))).^2)/params.sigma;
else
    diff_vec = 0;
    for c = 1:num_colors
        temp_I = I(:,:,c);
        diff_vec = diff_vec + ((temp_I(pairs(:,1))-temp_I(pairs(:,2))).^2)/params.sigma;
    end
    diff_vec = diff_vec/num_colors;
end

% form the weight matrix
num_pixels = nrows*ncols;
diff_vec = exp(-diff_vec); % ranges from 0 to 1 already

% fprintf('took [%3.3f seconds]\n\n',time_elapsed);
binary_cost_mat = sparse([pairs(:,1)],[pairs(:,2)],diff_vec,num_pixels,num_pixels);
time_elapsed = toc;
return;