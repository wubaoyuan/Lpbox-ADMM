function Reg = laicRegUpd(Pt, Reg0, diaSpeed)
% This function upadates each model point's trust region.
%
% Input
%   Pt        -  An Nm x 2 matrix. The ith row records ith model's
%                trust region center
%   Reg0      -  An Nm x 4 trust region matrix. The ith row records ith
%                model point's trust region. Each row has a data
%                structure like [x_min x_max y_min y_max]
%   diaSpeed  -  Trust region dismeter decreasing speed parameter
%
% Output
%   Reg       -  An Nm x 4 trust region matrix. The ith row records ith
%                model point's new trust region. Each row has a data
%                structure like [x_min x_max y_min y_max]
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 01-07-2012
%   modify    -  Feng Zhou (zhfe99@gmail.com), 07-28-2012

% dimension
Nm = size(Pt, 1);
Reg = zeros(size(Reg0));

% per point
for i = 1 : Nm
    x = Pt(i, 1);
    y = Pt(i, 2);
    
    xmin = Reg0(i, 1);
    xmax = Reg0(i, 2);
    ymin = Reg0(i, 3);
    ymax = Reg0(i, 4);

    diffXmin = abs(x - xmin);
    diffXmax = abs(xmax - x);

    diffYmin = abs(y - ymin);
    diffYmax = abs(ymax - y);
    
    % shrink in X
    if diffXmin + diffXmax > 2 * diaSpeed
        if diffXmin > diaSpeed && diffXmax > diaSpeed
            xmax = xmax - diaSpeed;
            xmin = xmin + diaSpeed;
        end
        if diffXmin < diaSpeed
            xmin = x + eps;
            xmax = xmax - (2 * diaSpeed - diffXmin);
        end
        if diffXmax < diaSpeed
            xmax = x - eps;
            xmin = xmin + (2 * diaSpeed - diffXmax);
        end
    end

    % shrink in Y
    if diffYmin + diffYmax > 2 * diaSpeed
        if diffYmin > diaSpeed && diffYmax > diaSpeed
            ymax = ymax - diaSpeed;
            ymin = ymin + diaSpeed;
        end
        if diffYmin < diaSpeed
            ymin = y + eps;
            ymax = ymax - (2 * diaSpeed - diffYmin);
        end
        if diffYmax < diaSpeed
            ymax = y - eps;
            ymin = ymin + (2 * diaSpeed - diffYmax);
        end
    end
    
    Reg(i, 1) = xmin;
    Reg(i, 2) = xmax;
    Reg(i, 3) = ymin;
    Reg(i, 4) = ymax;
end
