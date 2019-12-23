function ha = shGphCorLnC(gphs, X, parCor)
% Show correspondence as lines.
%
% Notice that the correspondence could contain continous values.
%
% Remark
%   n = n1 x n2
%
% Input
%   gphs    -  graphs, 1 x 2 (cell)
%     Pt    -  graph nodes, 2 x ni
%   X       -  correspondence, n1 x n2
%   parCor  -  parameter
%     th    -  threshold, {.3}
%     cls   -  color for line, {'k', 'b', 'g'} | 1 x 3 (cell)
%
% Output
%   ha
%     HCor  -  handle of correspondence line, k10 x k20 (cell)
%     cls   -  color for line, 1 x 2 (cell)
%     X     -  correspondence, k10 x k20
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-13-2012

% function option
th = ps(parCor, 'th', .3);
cls = ps(parCor, 'cls', {'k', 'b', 'g'});

% dimension
[n1, n2] = size(X);

% per pair of nodes
hold on;
HCor = cell(n1, n2);
for i = 1 : n1
    for j = 1 : n2
        % no correspondence
        if X(i, j) < th
            continue;
        end

        % coordinate
        lnX = [gphs{1}.Pt(1, i); gphs{2}.Pt(1, j)];
        lnY = [gphs{1}.Pt(2, i); gphs{2}.Pt(2, j)];

        % color and line width for incorrect correspondence
        lnWid = 1;
        cl = cls{1};

        % plot
        HCor{i, j} = plot(lnX, lnY, '-', 'Color', cl, 'LineWidth', lnWid);
    end
end

% store
ha.HCor = HCor;
ha.cls = cls;
ha.th = th;