
%% using KD-tree to compute the distance matrix
% run('D:\matlab code\vlfeat-0.9.19\toolbox\vl_setup')

t0 = cputime;
%% --------- Step 1 : generate the kd-tree
tic
kdtree = vl_kdtreebuild(dataset_matrix') ;
toc

%% -------- Step 2 : based on the kd-tree, compute the k-nearest neighbors and distances
num_batch = fix(num_sample/batch);
index_cell = cell(1,num_batch);
distance_cell = cell(1,num_batch);

if num_sample > 5000
        tic
        parfor iter_batch = 1:num_batch %1:batch:num_sample
            i = (iter_batch-1)*batch +1 : iter_batch*batch;
%              tic
                [index_cell{iter_batch}, distance_cell{iter_batch} ] = vl_kdtreequery(kdtree, dataset_matrix', dataset_matrix(i,:)', 'NumNeighbors', num_neighbor_size+1) ;    
 %             toc
        end
        toc
else
  %      tic
        for iter_batch = 1:num_batch %1:batch:num_sample
            i = (iter_batch-1)*batch +1 : iter_batch*batch;
            [index_cell{iter_batch}, distance_cell{iter_batch} ] = vl_kdtreequery(kdtree, dataset_matrix', dataset_matrix(i,:)', 'NumNeighbors', num_neighbor_size+1) ;    
        end
   %     toc
end

index_matrix = [ index_cell{:} ];
distance_matrix = [ distance_cell{:} ]; 

i = num_batch * batch + 1 : num_sample;
tic
[index_matrix_end, distance_matrix_end] = vl_kdtreequery(kdtree, dataset_matrix', dataset_matrix(i,:)', 'NumNeighbors', num_neighbor_size+1) ;    
toc
index_matrix = [ index_matrix, index_matrix_end];
distance_matrix = [ distance_matrix, distance_matrix_end];

index_matrix(1,:) = []; 
distance_matrix(1,:) = []; 
toc

%% normalize distance_matrix by dividing the maximal distance in the matrix, to (0, 1]
distance_max_value = max(distance_matrix(:));
distance_matrix_normalized = distance_matrix./ repmat(distance_max_value, num_neighbor_size, num_sample); 

%% square distance matrix 
distance_matrix_square = zeros(num_sample, num_sample);
for i = 1:num_sample
    distance_matrix_square( index_matrix(:,i) , i) = distance_matrix_normalized(:, i); 
end
distance_matrix_square = (distance_matrix_square + distance_matrix_square')./2;
