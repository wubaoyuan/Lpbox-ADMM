% Timothee Cour, 21-Apr-2008 17:31:23

% generate random sparse symmetric matrix
n = 1000;
W = sprandsym(n, .01);

% generate random affine constraint of the form Cx=b 
% (needs to be full rank with m < n)
m = 100;
C = rand(m, n);
b = rand(m, 1);

% compute to few eigevectors under affine constraint
nbEigenvectors = 10;
tic;
[X, lambda] = computeEigenvectorsAffineConstraint(W, C, b, nbEigenvectors);
toc;

% verify that constraint is satisfied (up to numerical precision, which can be set)
temp = C * X - repmat(b, 1, size(X,2));
constraintViolation = sqrt(sum(temp .^ 2, 1));

norm1(constraintViolation)
