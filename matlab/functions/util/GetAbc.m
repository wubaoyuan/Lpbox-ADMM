function [A,b,c] = GetAbc(unary, pairwise)
n = size(pairwise, 1);
b = unary(2,:)-unary(1,:); b=b(:);
const = sum(unary(1,:));
e=ones(n,1);
A  = -pairwise*4; 
We = -sum(A,2);

% compute the quadratic term efficiently
indx = sub2ind([n n],1:n,1:n);
A(indx) = A(indx)+We';

% get the constant
c=const;

return;
% fobj = x'*b + 0.5*x'*A*x + const;
