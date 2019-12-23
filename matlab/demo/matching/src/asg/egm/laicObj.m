function [obj, obj1, obj2, obj1s, obj2s] = laicObj(X, Pt2, C, H, lam)
% Compute the objective function of LAIC.
%
% Input
%   X       -  correspondence matrix, Nm x Nt
%   Pts     -  point position, 1 x 2 (cell), 2 x Nm, 2 x Nt
%   H       -  reconstruction matrix of the template point set, Nm x Nm
%   C       -  feature matching cost matrix, Nm x Nt
%
% Output
%   obj     -  objective
%   obj1    -  objective for 1st order term
%   obj2    -  objective for 2nd order term
%   obj1s   -  objective for 1st order term, 1 x Nm
%   obj2s   -  objective for 2nd order term, 1 x Nm
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-27-2012

if isempty(X)
    obj = 0;
    obj1 = 0;
    obj2 = 0;
    obj1s = [];
    obj2s = [];
    return;
end

% 1st order term
obj1s = diag(C * X');
obj1 = sum(obj1s);

% 2nd order term
obj2s = sum(abs(H * X * Pt2'), 2) * lam;
obj2 = sum(obj2s);

% obj
obj = obj1 + obj2;
