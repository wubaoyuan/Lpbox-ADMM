% function Y = multGXH(IndG, X, IndH)
%
% Compute the following matrix product in a much faster way.
%   G * X * H, where G and H are binary matrices,
%
% Input
%   IndG    -  index of non-zero variables in G, 2 x (lenG + 1)
%   X       -  nG x mH
%   IndH    -  index of non-zero variables in H, 2 x (lenH + 1)
%
% Output
%   Y       -  result, mG x nH
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

