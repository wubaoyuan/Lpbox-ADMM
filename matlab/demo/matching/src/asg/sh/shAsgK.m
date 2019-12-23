function shAsgK(K, KP, KQ, Ax)
% Show the correspondence matrix.
% 
% Input
%   K       -  affinity, n1n2 x n1n2 (sparse)
%   KP      -  node-node affinity, n1 x n2
%   KQ      -  edge-edge affinity, m1 x m2
%   Ax      -  axis, 1 x 3 (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-19-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-13-2011

if ~isempty(K)
    shM(K, 'ax', Ax{1});
    title('K');
end

shM(KP, 'ax', Ax{2});
title('KP');

shM(KQ, 'ax', Ax{3});
title('KQ');
