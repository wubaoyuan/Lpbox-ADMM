function shIts(objss, varargin)
% Show objective at each iteration.
%
% Input
%   objss    -  objective, 1 x m (cell), 1 x nIti
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
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% show option
psSh(varargin);

% function option
lnWid = ps(varargin, 'lnWid', 1);
mkSiz = ps(varargin, 'mkSiz', 0);
isMkEg = psY(varargin, 'mkEg', 'n');
algs = ps(varargin, 'algs', []);

% dimension
m = length(objss);
nIts = cellDim(objss, 2);

hold on;

% plot obj
for j = 1 : m
    [~, cl] = genMkCl(j);
    plot(1 : nIts(j), objss{j}, '-', 'LineWidth', lnWid, 'Color', cl);
end

% legend
if ~isempty(algs)
    legend(algs{:});
end

% axis
axis square;
objs = mcat('horz', objss);
mi = min(objs);
ma = max(objs);;
ylim([mi, ma + (ma - mi) * .1]);
