function h = show3dUpd2(h0, vis, varargin)
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
c = ps(varargin, 'c', []);
lnWid = ps(varargin, 'lnWid', 1);
mkSiz0 = ps(varargin, 'mkSiz0', 2);
mkSiz1 = ps(varargin, 'mkSiz1', 5);
isFace = psY(varargin, 'face', 'n');
n = length(vis0);

h = h0;
h.vis = vis;

% check whether vis == vis0
if all(vis == vis0)
    return;
end

% mks
[mks, cls] = genMkCl;

% default color
cl0 = [1 1 1] * .5;
mk0 = mks{1};
 
% focused string in the current step
for i = 1 : n
    
    if vis(i)
        if isempty(c)
            cl = cls{i};
            mk = mks{i};
        else
            cl = cls{c};
            mk = mks{c};
        end
        mkSiz = mkSiz1;
    else
        cl = cl0;
        mk = mk0;
        mkSiz = mkSiz0;
        continue;
    end
    
    set(h0.poi{i}, 'Marker', mk, ...
        'MarkerSize', mkSiz, 'Color', cl, 'LineWidth', lnWid);

    if isFace
        set(h0.poi{i}, 'MarkerFaceColor', 'none');
    end
end
