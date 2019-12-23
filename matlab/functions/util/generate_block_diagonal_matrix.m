%% makes a block diagonal matrix from one matrix
% A is a square matrix

function [A_mat] = generate_block_diagonal_matrix(A,num_blocks)
[m n] = size(A);
A_mat = sparse(m*num_blocks,n*num_blocks);

for b = 1:num_blocks
    A_mat(1+(b-1)*m:b*m,1+(b-1)*n:b*n) = A;
end
return;