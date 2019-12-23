function [A,b,c] = convert2Abc(unary,pairwise)
n = size(pairwise, 1);

% get the unary term
b = unary(2,:)-unary(1,:); b=b(:);
const = sum(unary(1,:));

% get the binary term
A  = -pairwise;
We = -sum(A,2);

% compute the quadratic term efficiently
% ++ assume that the diagonal of the pariwise matrix is ALL zeros
A = A+sparse([1:n],[1:n],We',n,n); % add to the diagonal
A = 2*A;


% indx = sub2ind([n n],1:n,1:n);
% A(indx) = A(indx)+We';

% get the constant
c=const;
return;
