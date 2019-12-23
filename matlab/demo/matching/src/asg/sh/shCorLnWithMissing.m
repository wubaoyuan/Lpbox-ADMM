function ha = shCorLnWithMissing(Pts, P, PT, varargin)
% Show correspondence as lines.
%
% Remark
%   k0 = k01 x k02
%   If visPts{1} == [], then visPts{1} is default to ones(ki, 1)
%   If visPts{1} ~= [], then ki = length(find(visPts{i}))
%
% Input
%   Pts       -  graph node, 1 x 2 (cell), 2 x ki
%   P         -  correspondence, k10 x k20
%   PT        -  ground-truth correspondence, [] | k10 x k20
%   varargin
%     cls     -  color for line, {'k', 'b'} | 1 x 2 (cell)
%     visPts  -  node existence status, {[]} | 1 x 2 (cell), ki0 x 1
%
% Output
%   ha
%     HCor    -  handle of correspondence line, k10 x k20 (cell)
%     cls     -  color for line, 1 x 2 (cell)
%     P       -  correspondence, k10 x k20
%     PT      -  ground-truth correspondence, k10 x k20
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 02-06-2012

% function option
cls = ps(varargin, 'cls', {'k', 'b'});
visPts = ps(varargin, 'visPts', []);

% dimension
k1 = size(Pts{1}, 2);
k2 = size(Pts{2}, 2);
[k10a, k20a] = size(P);

% default visPts
if isempty(visPts)
    k10 = k1;
    k20 = k2;
    visPt1 = ones(k1, 1);
    visPt2 = ones(k2, 1);

else
    visPt1 = visPts{1};
    visPt2 = visPts{2};
    k10 = size(visPt1, 1);
    k20 = size(visPt2, 1);
end
% equal('k0s', [k10 k20], [k10a k20a]);

% index
idxPt1 = find(visPt1);
ordPt1 = zeros(1, k10);
ordPt1(idxPt1) = 1 : length(idxPt1);
idxPt2 = find(visPt2);
ordPt2 = zeros(1, k20);
ordPt2(idxPt2) = 1 : length(idxPt2);

% default ground-truth
if isempty(PT)
    PT = zeros(k10, k20);
end

% per pair
HCor = cell(k10, k20);
for i = 1 : k10
    for j = 1 : k20
        % no correspondence or missing point
        if visPt1(i) == 0 || visPt2(j) == 0 || (P(i, j) == 0 && PT(i, j) == 0)
            continue;
        end

        % coordinate
        lnX = [Pts{1}(1, ordPt1(i)); Pts{2}(1, ordPt2(j))];
        lnY = [Pts{1}(2, ordPt1(i)); Pts{2}(2, ordPt2(j))];
        %lnX = [Pts{1}(1, i); Pts{2}(1, j)];
        %lnY = [Pts{1}(2, i); Pts{2}(2, j)];

        % color
        % wrong correspondence
        if P(i, j) == 0 && PT(i, j) == 1
            lnWid = 2;
            cl = cls{2};
            
        % correct correspondence
        else
            lnWid = 1;
            cl = cls{1};
        end

        % plot
        HCor{i, j} = plot(lnX, lnY, '-', 'Color', cl, 'LineWidth', lnWid);
    end
end

% store
ha.HCor = HCor;
ha.cls = cls;
ha.P = P;
ha.PT = PT;
