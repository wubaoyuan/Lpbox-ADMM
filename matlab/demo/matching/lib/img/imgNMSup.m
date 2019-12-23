function [R, idx] = imgNMSup(R0, siz)
% Non-maximum suppression.
%
% Input
%   R0      -  original image, h x w
%   siz     -  size, 1 x 1
%
% Output
%   R       -  new image, h x w
%   idx     -  index of the peak point, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-10-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-23-2012

% dimension
[h, w] = size(R0);
n = h * w;
maskH = siz;
maskW = siz;
m = maskH * maskW;
 
[I, J] = ndgrid(1 : h, 1 : w);
I = repmat(I(:)', [m 1]);
J = repmat(J(:)', [m 1]);
 
mh = (maskH - 1) / 2;
mw = (maskW - 1) / 2;
[K, L] = ndgrid(-mh : mh, -mw : mw);
K = repmat(K(:), [1 n]);
L = repmat(L(:), [1 n]);
 
I = I + K; I = I(:);
J = J + L; J = J(:);
index = sub2ind_modify([h w], I, J);
Index = reshape(index, m, n);
 
% get the max for each window
R1 = [R0(:); eps];
im_max = R1(index);
im_max = reshape(im_max, m, n);
[im_max, p] = max(im_max, [], 1);
ind = sub2ind([m, n], p, 1 : n);
idx = index(ind)';
idx0 = 1 : n;
 
% suppress non-local maxima to zero
R = R0;
vis = idx0 == idx;
R(~vis) = 0;
idx = find(idx0 == idx);

% sort
vs = R(idx);
[~, ord] = sort(vs, 'descend');
idx = idx(ord);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = sub2ind_modify(matrixSize, rowSub, colSub)
index = ones(size(rowSub)) * (matrixSize(1) * matrixSize(2) + 1);
 
idx = rowSub >= 1 & rowSub <= matrixSize(1) ...
      & colSub >= 1 & colSub <= matrixSize(2);
index(idx) = sub2ind(matrixSize, rowSub(idx), colSub(idx));
