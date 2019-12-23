function setAx(box, parAx)
% Set properties of current axes.
%
% Input
%   box      -  bounding box, d x 2
%   parAx    -  parameter
%     set    -  setting flag, {'y'} | 'n'
%     grid   -  grid on flag, 'y' | {'n'}
%     eq     -  axis equal flag, 'y' | {'n'}
%     sq     -  axis square flag, 'y' | {'n'}
%     ij     -  axis ij flag, 'y' | {'n'}
%     ax     -  axis flag, {'y'} | 'n'
%     tick   -  tick flag, {'y'} | 'n'
%     label  -  label flag, 'y' | {'n'}
%     ang    -  view angle, {[]}
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 01-29-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function parameter
isSet = psY(parAx, 'set', 'y');
isGrid = psY(parAx, 'grid', 'n');
isEq = psY(parAx, 'eq', 'n');
isSq = psY(parAx, 'sq', 'n');
isIj = psY(parAx, 'ij', 'n');
isAx = psY(parAx, 'ax', 'y');
isTick = psY(parAx, 'tick', 'y');
isLabel = psY(parAx, 'label', 'n');
ang = ps(parAx, 'ang', [30 80]);

% bounding box
d = size(box, 1);

% adjust axis
if isSet
    % grid
    if isGrid
        grid on;
    else
        grid off;
    end

    % axis
    if isEq
        axis equal;
    end
    if isSq
        axis square;
    end
    if isIj
        axis ij;
    end
    if isAx
        axis on;
    else
        axis off;
    end

    % limitation
    if abs(box(1, 1) - box(1, 2)) > eps
        xlim([box(1, 1), box(1, 2)]);
    end
    if abs(box(2, 1) - box(2, 2)) > eps
        ylim([box(2, 1), box(2, 2)]);
    end
    if d == 3 && abs(box(3, 1) - box(3, 2)) > eps
        zlim([box(3, 1), box(3, 2)]);
    end
    
    % tick
    if ~isTick
        if d == 3
            set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
        else
            set(gca, 'XTick', [], 'YTick', []);
        end
    end
    
    % label
    if isLabel
        if d == 3
            xlabel('x');
            ylabel('y');
            zlabel('z');
        else
            xlabel('x');
            ylabel('y');
        end
    end

    % view angle
    if d == 3 && ~isempty(ang)
        view(ang);
    end
end
