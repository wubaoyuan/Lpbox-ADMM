function Reg = laicRegSet(Pt, dia)
% Sets initial trust regions for each model points. 
% 
% Input
%   Pt      -  The ith row records ith model's trust region center, n x 2
%   dia     -  The initial trust region diameter
%           
% Output    
%   Reg     -  trust region matrix, n x 4
%              Each row: [x_min x_max y_min y_max]
%           
% History   
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-16-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-07-2012

% dimension
n = size(Pt, 1);

% per point
Reg = size(n, 4);
for i = 1 : n
    x = Pt(i, 1);
    y = Pt(i, 2);
    
    Reg(i, 1) = x - dia;
    Reg(i, 2) = x + dia;
    Reg(i, 3) = y - dia;
    Reg(i, 4) = y + dia;
end
