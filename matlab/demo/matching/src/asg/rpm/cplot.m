function cplot(in1, in2, in3, in4, in5, in6)
% Plot points.
%
% Usage: 
% [] = cplot (x);
% [] = cplot (x,y);
% [] = cplot (x,y,z);
%
% [] = cplot (x, marker_str, marker_size);
% [] = cplot (x, xmarker_str, xmarker_size, y, ymarker_str, y_marker_size);
%
% Input
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012

% check input:
if (nargin == 1) % -------------------------------- (x) ---
    x = in1; xmarker_str = 'go'; xmarker_size = 6;
    y = [];  ymarker_str = 'r+'; ymarker_size = 6;
    z = [];  zmarker_str = 'bo'; zmarker_size = 6;
    
elseif (nargin == 2) % -------------------------- (x,y) ---
    x = in1; xmarker_str = 'go'; xmarker_size = 6;
    y = in2; ymarker_str = 'r+'; ymarker_size = 6;
    z = [];  zmarker_str = 'bo'; zmarker_size = 6;
    
elseif (nargin == 3) & (~isstr(in2)) % -------- (x,y,z) ---
    x = in1; xmarker_str = 'go'; xmarker_size = 6;
    y = in2; ymarker_str = 'r+'; ymarker_size = 6;
    z = in3; zmarker_str = 'bo'; zmarker_size = 6;
    
elseif (nargin == 3) & (isstr(in2)) % ---- (x, 'go', 3) ---
    x = in1; xmarker_str = in2; xmarker_size = in3;
    y = [];
    z = [];
    
elseif (nargin == 6) % -------- (x, 'go, 3, y, 'r+', 3) ---
    x = in1; xmarker_str = in2; xmarker_size = in3;
    y = in4; ymarker_str = in5; ymarker_size = in6;
    z = [];
    
else
    error('# ERROR #: cplot -- wrong input!');
end

% plot x:
[n, dim] = size(x);
if n >= 1
    cplot_1pointset(x, xmarker_str, xmarker_size); 
end
hold on;

% plot y:
[n, dim] = size(y);
if n >= 1
    cplot_1pointset(y, ymarker_str, ymarker_size); 
end

% plot z:
[n, dim] = size(z);
if n >= 1
    cplot_1pointset(z, zmarker_str, zmarker_size); 
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cplot_1pointset(x, xmarker_str, xmarker_size)
% Plot one point set.

[n, dim] = size(x);
if (dim == 2)
    h = plot(x(:,1), x(:,2), xmarker_str, 'markersize', xmarker_size); 
    axis('equal');
    hold on;
    
elseif (dim == 3)
    h = plot3(x(:,1), x(:,2), x(:,3), xmarker_str, 'markersize', xmarker_size); 
    axis('equal'); 
    set(gca, 'box', 'on');
    hold on;
end;
