function rpmPlot2(Ax, x, y, vx, m, thM, T, algT, c_tps, d_tps, w, sigma, algX)
% Plot the cMIX progress. (no Cx, just T -- as isotropic covariance).
%
% Input
%   Ax     -  Axes
%   x      -  1st point set
%   y      -  2nd point set
%   vx
%   m
%   thM
%   T
%   algT
%   c_tps
%   d_tps
%   w
%   sigma
%   algX
%
% History
%   create  -  Anand Rangarajan (anand@noodle.med.yale.edu), 04-27-2000
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-23-2012

% configure
xmarker = 'go'; 
xsize = 6;
ymarker = 'r+'; 
ysize = 3;

% dimension
[siz1, ~] = size(x);
[siz2, dim] = size(y);
c = c_tps;
d = d_tps;

% original point
set(gcf, 'CurrentAxes', Ax{1}); 
cla;
cplot(x, xmarker, xsize);
hold on; 
cplot(y, ymarker, ysize); 
axis('equal');
title('Original V and X'); 

% transformed point
set(gcf, 'CurrentAxes', Ax{2});
cla;
cplot(vx, xmarker, xsize);
hold on;
cplot(y, ymarker, ysize); 
axis('equal');
title('Transformed V + X');

% circle
set(gcf, 'CurrentAxes', Ax{3});
cla;
cplotg(vx, y, m, thM);
hold on;
if ~strcmp(algX, 'icp')
    cMIX_plot_mixture_simple(vx, T, 'b-'); 
    hold on;
end
cplot(vx, xmarker, xsize);
hold on;
cplot(y, ymarker, ysize);
hold on;
axis('equal');
title('Transformed V + X');

% TPS
set(gcf, 'CurrentAxes', Ax{4});
cla;
switch(algT)
  case 'tps'
    ctps_plot_grid(x, x, c, d);
    axis('equal');
    title ('TPS Warping');

  case 'rbf'
    crbf_plot_grid(x, x, w, sigma); 
    axis('equal'); 
    title('RBF Warping');
end
cplot(y, ymarker, ysize); 
hold on;
cplot(vx, xmarker, xsize);
hold on;
axis('equal');

% estimated
set(gcf, 'CurrentAxes', Ax{5});
cla;
vy  = m * y ./  ((sum(m'))' * ones(1,dim));
cplot(y, 'r.', ysize); 
hold on;
cplot(vy, ymarker, xsize); 
axis('equal');  
title('Estimated Shape Y=MX');

shM(m, 'ax', Ax{6});
title('correspondence');