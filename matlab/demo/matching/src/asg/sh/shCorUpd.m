function ha = shCorUpd(ha, Pts, varargin)
% Update correspondence.
%
% Input
%   ha        -  handle
%   Pts       -  graph node, 1 x 2 (cell), 2 x ki
%   varargin
%     visPts  -  node existence status, {[]} | 1 x 2 (cell), ki0 x 1
%
% Output
%   ha        -  handle
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
visPts = ps(varargin, 'visPts', []);

% dimension
ks = cellDim(Pts, 2);

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

% index
idxPt1 = find(visPt1);
ordPt1 = zeros(1, k10);
ordPt1(idxPt1) = 1 : length(idxPt1);
idxPt2 = find(visPt2);
ordPt2 = zeros(1, k20);
ordPt2(idxPt2) = 1 : length(idxPt2);

% per pair
for i = 1 : ks(1)
    for j = 1 : ks(2)
        % no correspondence or missing data
        if ha.P(i, j) == 0 || visPt1(i) == 0 || visPt2(j) == 0
            if ~isempty(ha.HCor{i, j})
                delete(ha.HCor{i, j});
                ha.HCor{i, j} = [];
            end
            continue;
        end
        
        % coordinate
        lnX = [Pts{1}(1, ordPt1(i)); Pts{2}(1, ordPt2(j))];
        lnY = [Pts{1}(2, ordPt1(i)); Pts{2}(2, ordPt2(j))];

        % color
        if ha.P(i, j) == ha.PT(i, j)
            lnWid = 1;
            cl = ha.cls{1};
        else
            lnWid = 2;
            cl = ha.cls{2};
        end
        
        % plot
        if isempty(ha.HCor{i, j})
            ha.HCor{i, j} = plot(lnX, lnY, '-', 'Color', cl, 'LineWidth', lnWid);        
        else
            set(ha.HCor{i, j}, 'XData', lnX, 'YData', lnY);
        end
    end
end
