function D = conDst(X1, X2, varargin)
% Compute (squared) distance matrix.
%
% Remark
%   Dij is the squared Euclidean distance between the i-th point in X1 and j-th point in X2.
%   i.e., D(i,j) = || X1(:, i) - X2(:, j) ||_2^2
%
% Usage (1)
%   input   -  X1 = rand(3, 5); X2 = rand(3, 6);
%   call    -  D = conDst(X1, X2);
%
% Usage (2): saving the D in the global workspace
%   input   -  X1 = rand(3, 5); X2 = rand(3, 6);
%   call    -  conDst(X1, X2);
%
% Input
%   X1      -  1st sample matrix, dim x n1
%   X2      -  2nd sample matrix, dim x n2
%   varargin
%     dst   -  distance type, {'e'} | 'b'
%              'e': Euclidean distance
%              'b': binary distance
%
% Output
%   D       -  squared distance matrix, n1 x n2
%              If nargout == 0, then the distance matrix is saved in the global variable (DG).
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-25-2012

% global
global DG;
isDG = nargout == 0; 

% function option
dst = ps(varargin, 'dst', 'e');

% dimension
n1 = size(X1, 2);
n2 = size(X2, 2);
if size(X1, 1) == 1
    X1 = [X1; zeros(1, n1)]; 
    X2 = [X2; zeros(1, n2)]; 
end
XX1 = sum(X1 .* X1); XX2 = sum(X2 .* X2); 

% compute
if isDG
    DG = -2 * X1' * X2;
    for i2 = 1 : n2
        DG(:, i2) = DG(:, i2) + XX1';
    end
    for i1 = 1 : n1
        DG(i1, :) = DG(i1, :) + XX2;
    end
else
    X12 = X1' * X2;
    D = repmat(XX1', [1, n2]) + repmat(XX2, [n1, 1]) - 2 * X12;
end

% Euclidean distance
if strcmp(dst, 'e')

% binary distance
elseif strcmp(dst, 'b')
    if isDG
        DG = real(DG > 1e-8);
    else
        D = real(D > 1e-8);
    end

else
    error(['unknown distance type: ' dst]);
end
