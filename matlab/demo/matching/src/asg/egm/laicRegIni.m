function Reg = laicRegIni(X, Pt, bas)
% Sets initial trust regions for each model points. 
% 
% Input
%   X       -  correspondence matrix, k x n
%   Pt      -  2 x n
%   bas     -  basis, 1 x k (cell), 1 x ni
%           
% Output    
%   Reg     -  region matrix, 2 x 2 x k
%           
% History   
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-16-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-28-2012

% dimension
[k, n] = size(X);

% center 
Cen = Pt * X';

% per point
Reg = zeros(2, 2, k);
for c = 1 : k
    cen = Cen(:, c);
    
    Ptc = Pt(:, bas{c});

    boxMi = min(Ptc, [], 2);
    boxMa = max(Ptc, [], 2);
    
    siz1 = cen - boxMi;
    siz2 = boxMa - cen;
    siz = max([siz1, siz2], [], 2);
    
    Reg(:, :, c) = cen(:, [1 1]) + [-siz, siz];
end
