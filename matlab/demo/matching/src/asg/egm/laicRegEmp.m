function vis = laicRegEmp(Reg, T)
% Tests if any trust region is empty and outputs empty indices.
%
% Input
%   Reg     -  Nm x 4 trust region matrix. The ith row records ith
%              model point's trust region. Each row has a data
%              structure like [x_min x_max y_min y_max].
%   T       -  An Nt x 2 matrix recording 2D Nt target points' 
%              coordinates.
%           
% Output    
%   vis     -  An Nm x 1 0-1 vector. The ith element represents
%              whether the ith trust region has no point in it.
%           
% History   
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-07-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-07-2012

% dimension
Nm = size(Reg, 1);
Nt = size(T, 1);
vis = false(Nm, 1);

for i = 1 : Nm
    xmin = Reg(i, 1);
    xmax = Reg(i, 2);
    ymin = Reg(i, 3);
    ymax = Reg(i, 4); 
    
    xminAll = repmat(xmin, Nt, 1);
    xmaxAll = repmat(xmax, Nt, 1);
    yminAll = repmat(ymin, Nt, 1);
    ymaxAll = repmat(ymax, Nt, 1);
    
    isOutAll = T(:, 1) < xminAll | T(:, 1) > xmaxAll | T(:, 2) < yminAll | T(:, 2) > ymaxAll;
    vis(i) = all(isOutAll);
end
