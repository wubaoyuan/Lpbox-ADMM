function h = show3dUpd(h0, vis, varargin)
% Show a string list (Updating).
%
% Input
%   h0      -  original handle
%   vis     -  position of new focused string, 1 x n
%   varargin
%     cl    -  string color, {[1 1 0]}
%
% Output
%   h       -  new handle
%     vis   -  position of new focused string, 1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-31-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% focused string in the previous step
vis0 = ps(h0, 'vis', []);
n = length(vis0);

h = h0;
h.vis = vis;

% check whether vis == vis0
if all(vis == vis0)
    return;
end
 
for i = 1 : n
    if vis0(i)
        set(h0.poi{i}, 'visible', 'on');
    else
        set(h0.poi{i}, 'visible', 'off');
    end
end

% focused string in the current step
for i = 1 : n
    if vis(i)
        set(h0.poi{i}, 'visible', 'on');
    else
        set(h0.poi{i}, 'visible', 'off');
    end
end
