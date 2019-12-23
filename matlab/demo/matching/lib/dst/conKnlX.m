function S = conKnlX(X, para)
% Compute kernel matrix from sample matrix
%
% Input
%   X       -  sample matrix, dim x n
%   para    -  parameter
%     dst   -  see function conDist for details
%     knl   -  see function conSim for details
%     nei   -  see function conSim for details
%
% Output
%   S       -   similarity matrix, n x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

D = conDist(X, X, 'dst', para.dst);
S = conSim(D, 'knl', para.knl, 'nei', para.nei);
