function ha = shGphCorLn(gphs, X, XT, parCor)
% Show correspondence as lines.
%
% Remark
%   n = n1 x n2
%
% Input
%   gphs    -  graphs, 1 x 2 (cell)
%     Pt    -  graph nodes, 2 x ni
%   X       -  correspondence, n1 x n2
%   XT      -  ground-truth correspondence, [] | n1 x n2
%   parCor  -  parameter
%     cls   -  color for line, {'k', 'b', 'g'} | 1 x 3 (cell)
%
% Output
%   ha
%     HCor  -  handle of correspondence line, k10 x k20 (cell)
%     cls   -  color for line, 1 x 2 (cell)
%     X     -  correspondence, k10 x k20
%     XT    -  ground-truth correspondence, k10 x k20
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 11-20-2012

% function option
cls = ps(parCor, 'cls', {'k', 'b', 'g'});
cls = ps(parCor, 'cls', {'g', 'b', 'b'});
lnWidIn = 1;
lnWidC = 1;

% dimension
[n1, n2] = size(X);

% default ground-truth
if isempty(XT)
    XT = zeros(n1, n2);
end

% per pair
HCor = cell(n1, n2);
for i = 1 : n1
    for j = 1 : n2
        % no correspondence
        if X(i, j) == 0 && XT(i, j) == 0
            continue;
        end

        % coordinate
        lnX = [gphs{1}.Pt(1, i); gphs{2}.Pt(1, j)];
        lnY = [gphs{1}.Pt(2, i); gphs{2}.Pt(2, j)];

        % color and line width for incorrect correspondence
        if X(i, j) == 0 && XT(i, j) == 1
            lnWid = lnWidIn;
            cl = cls{2};

        elseif X(i, j) == 1 && XT(i, j) == 0
            lnWid = lnWidC;
            cl = cls{3};

        % color and line width for correct correspondence
        else
            lnWid = lnWidC;
            cl = cls{1};
        end

        % plot
        HCor{i, j} = plot(lnX, lnY, '-', 'Color', cl, 'LineWidth', lnWid);
    end
end

% store
ha.HCor = HCor;
ha.cls = cls;
ha.X = X;
ha.XT = XT;
