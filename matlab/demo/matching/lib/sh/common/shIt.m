function shIt(Obj, its, varargin)
% Show objective at each iteration.
%
% Input
%   Obj      -  objective, m x nIt
%   its      -  iteration ids, 1 x nIt
%   varargin
%     show option
%     axis   -  axis flag, {'y'} | 'n'
%     lnWid  -  line width, {1}
%     mkSiz  -  size of mks, {0}
%     mkEg   -  flag of marker edge, 'y' | {'n'}
%     lim0   -  predefined limitation, {[]}
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 09-20-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 11-08-2012

% show option
psSh(varargin);

% function option
lnWid = ps(varargin, 'lnWid', 1);
mkSiz = ps(varargin, 'mkSiz', 0);
isMkEg = psY(varargin, 'mkEg', 'n');
itNms = ps(varargin, 'itNms', []);
itMa = ps(varargin, 'itMa', 2);

% dimension
[m, nIt] = size(Obj);

% iterations
if isempty(its)
    its = ones(1, nIt);
end
if isempty(itNms)
    k = max(its);
else
    k = length(itNms);
end
vis = zeros(1, k);

hold on;

% per iteration
for c = 1 : k
    idx = find(its == c);
    if isempty(idx)
        continue;
    end
    vis(c) = 1;

    [mk, cl] = genMkCl(c);
    for j = 1 : m
        hTmp = plot(idx, Obj(j, idx), mk, 'Color', cl);

        if mkSiz > 0
            set(hTmp, 'Marker', mk, 'MarkerSize', mkSiz, 'MarkerFaceColor', cl);
            if isMkEg
                set(hTmp, 'MarkerEdgeColor', 'k');
            end
        end
    end
end
if ~isempty(itNms)
    legend(itNms(vis == 1));
end

% plot the line
for j = 1 : m
    [~, cl] = genMkCl(j);
    plot(1 : nIt, Obj(j, :), '-', 'LineWidth', lnWid, 'Color', cl);
end

% axis
axis square;
if itMa == 0
    ma = max(Obj(:));
elseif nIt >= itMa
    ma = max(Obj(:, itMa));
end
mi = min(Obj(:));
if nIt > itMa && ma > eps
    if ma > mi + eps
        ylim([mi, ma + (ma - mi) * 1.1]);
    else
        ylim([mi, ma * 1.1]);
    end
end
