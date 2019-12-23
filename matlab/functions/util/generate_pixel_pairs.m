%% finds the linear indices of pixels neighboring each pixel and saves them as pairs
% + k defines a (2k+1)x(2k+1) centered at the pixel

function [pairs,offsets,neighborhood_mat] = generate_pixel_pairs(nrows,ncols,k)
% get offsets in the spatial indices
[offsets] = get_offsets(k);

% get all spatial indices
[rmat cmat] = meshgrid(1:nrows,1:ncols);
spatial_indices = [rmat(:) cmat(:)];

num_pixels = numel(rmat);
neighborhood_mat = [];
for p = 1:size(offsets,1)
    offset = offsets(p,:);
    if norm(offset)>0
        neighborhood_mat = [neighborhood_mat; [spatial_indices spatial_indices+repmat(offset,num_pixels,1)]];
    end
end


% keep only valid indices
valid_idx = (neighborhood_mat(:,3)>=1) & (neighborhood_mat(:,3)<=nrows) & (neighborhood_mat(:,4)>=1) & (neighborhood_mat(:,4)<=ncols);
pairs = [sub2ind([nrows ncols],neighborhood_mat(valid_idx,1),neighborhood_mat(valid_idx,2)) sub2ind([nrows ncols],neighborhood_mat(valid_idx,3),neighborhood_mat(valid_idx,4))];
neighborhood_mat = neighborhood_mat(valid_idx,:);
return;

%% gets the linear indices of a point and its neighbors defined by k
function [offsets] = get_offsets(k)

row_offsets = repmat([-k:1:k]',1,2*k+1);
col_offsets = row_offsets';
offsets = [row_offsets(:) col_offsets(:)];

return

%% vectorize matrix
function [a] = vec(A)
a = A(:);
return;

