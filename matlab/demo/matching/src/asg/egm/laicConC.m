function C = laicConC(Rs, Pt)
% Compute the cost matrix for LAIC algorithm.
%
% Input
%   Rs      -  response, h x w x k
%   Pt      -  testing point set, 2 x n
%
% Output
%   C       -  cost matrix, k x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-20-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-26-2012

% dimension
k = size(Rs, 3);
n = size(Pt, 2);

% node similarity
K = zeros(k, n);
for iTr = 1 : k
    Ri = Rs(:, :, iTr);
    Ri(Ri(:) < 0) = 0; % maybe a bug, not a good fix
    
    K(iTr, :) = log(imgPtVal(Ri, Pt)');
end

% simiarlity -> distance
C = -K;

% remove inf
idx = find(C == inf);
C(idx) = 1e8;
