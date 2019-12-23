function [D, DSum] = conDstMask(Pts, varargin)
% Construct the distane of binary masks.
%
% Input
%   Pts       -  binary masks, 1 x nF (cell)
%
% Output
%   D         -  distance matrix, nF x nF
%   DSum      -  area size, nF x nF
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
nF = length(Pts);
Lns = cell(1, nF);
for iF = 1 : nF
    M = maskP2M([], Pts{iF});
    Lns{iF} = maskM2L(M);
end

% area size
Area = zeros(1, nF);
for iF = 1 : nF
    Ln = Lns{iF};
    Area(iF) = sum(Ln(3, :));
end
DSum = zeros(nF, nF);
for i = 1 : nF
    for j = 1 : nF
        DSum(i, j) = Area(i) + Area(j);
    end
end

% compare
D = zeros(nF, nF);
for i = 1 : nF
    for j = i + 1 : nF
        D(i, j) = 1 - maskMatchSearch(Lns{j}, Lns{i}, [.2; .2]);
    end
end
D = D + D';
