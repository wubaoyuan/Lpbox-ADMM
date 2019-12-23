function setAxTick(pos, form, ticks, vals)
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

% dimension
n = length(ticks);

labels = cell(1, n);
for i = 1 : n
    labels{i} = sprintf(form, vals(i));
end

if strcmp(pos, 'x')
    set(gca, 'XTick', ticks, 'XTickLabel', labels);
elseif strcmp(pos, 'y')
    set(gca, 'YTick', ticks, 'YTickLabel', labels);
else
    error('unknown position: %s', pos);
end

    
