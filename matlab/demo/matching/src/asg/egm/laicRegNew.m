function [bas, Reg] = laicRegNew(X, Pt, ba0s, Reg0, rat)
% Sets initial trust regions for each model points. 
% 
% Input
%   X       -  correspondence, k x n
%   Pt      -  2 x n
%   ba0s    -  original basis, 1 x k (cell), 1 x ni
%   Reg0    -  original region matrix, 2 x 2 x k
%   rat     -  ratio
%           
% Output    
%   bas     -  basis, 1 x k (cell), 1 x ni
%   Reg     -  region matrix, 2 x 2 x k
%           
% History   
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-16-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-16-2012

% dimension
k = length(ba0s);

% center 
Cen = Pt * X';

% per point
Reg = zeros(2, 2, k);
bas = cell(1, k);
for c = 1 : k
    Reg0c = Reg0(:, :, c);
    
    % original size
    siz0 = Reg0c(:, 2) - Reg0c(:, 1);
    
    % reduce the size
    siz = siz0 .* rat;
    
    % new center
    cen = Cen(:, c);
    
    % new region
    Regc = cen(:, [1 1]) + [-siz, siz] / 2;
    
    % original basis
    Pt0c = Pt(:, ba0s{c});
    
    % inside the region
    vis = Pt0c(1, :) >= Regc(1, 1) & Pt0c(2, :) >= Regc(2, 1) & ...
          Pt0c(1, :) <= Regc(1, 2) & Pt0c(2, :) <= Regc(2, 2);
    
    % store
    Reg(:, :, c) = Regc;
    bas{c} = ba0s{c}(vis);
end
