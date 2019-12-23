function rpmPlot(Ax, method, x, y, z, vx, m, m_threshold, T, algT, c_tps, d_tps, w, sigma_kernel, m_method)
% Plot the cMIX progress. (no Cx, just T -- as isotropic covariance).
%
% Input
%   Ax
%   method, 2 | ...
%   x
%   y
%   z
%   vx
%   m
%   m_threshold
%   T
%   algT
%   c_tps
%   d_tps
%   w
%   sigma_kernel
%   m_method
    
% configure
xmarker = 'go'; xsize = 6;
ymarker = 'r+'; ysize = 3;
zmarker = 'go'; zsize = 6;
inter_marker = 'g+'; 

% dimension
[siz1, dim] = size(x);
[siz2, dim] = size(y);
c = c_tps;
d = d_tps;

if (dim == 2)
  
switch(method)
  case 2
    set(gcf, 'CurrentAxes', Ax{1}); 
    cla;
    cplot(x, xmarker, xsize);
    hold on; 
    cplot(y, ymarker, ysize); 
    axis('equal');
    title('Original V and X'); 
    
    set(gcf, 'CurrentAxes', Ax{2});
    cla;
    cplot(vx, xmarker, xsize); 
    hold on; 
    cplot(y, ymarker, ysize); 
    axis('equal');
    title('Transformed V + X');

    set(gcf, 'CurrentAxes', Ax{3});
    cla;
    cplotg(vx, y, m, m_threshold); 
    hold on;
    if ~strcmp(m_method, 'icp')
        cMIX_plot_mixture_simple(vx, T, 'b-'); 
        hold on;
    end
    cplot(vx, xmarker, xsize);
    hold on; 
    cplot(y, ymarker, ysize);
    hold on;
    axis('equal');
    title('Transformed V + X');

    set(gcf, 'CurrentAxes', Ax{4});
    cla;
    switch(algT)
      case 'tps'
	ctps_plot_grid(x, x, c, d);
	axis('equal');
        title ('TPS Warping');
        
      case 'rbf'
	crbf_plot_grid(x, z, w, sigma_kernel); 
	axis('equal'); 
        title('RBF Warping');
    end
    cplot(y, ymarker, ysize); 
    hold on;
    cplot(vx, xmarker, xsize);
    hold on;
    axis('equal');
	
    set(gcf, 'CurrentAxes', Ax{5});
    cla;
    vy  = m * y ./  ((sum(m'))' * ones(1,dim));
    cplot (y, 'r.', ysize); hold on;
    cplot (vy, ymarker, xsize); 
    axis('equal');  
    title ('Estimated Shape Y=MX');
  
  case 0
    cplotg (vx, y, m, m_threshold); hold on; 
    cplot  (vx, xmarker, xsize); hold on; 
    cplot  (y,  ymarker, ysize); hold on; 
    cMIX_plot_mixture_simple (vx, T); title ('Transformed V + X');
    
  case 1
    
    h_sub1 = subplot ('position', [0.05 0.6 0.25 0.3]); 
    cplot (x, xmarker, xsize); hold on; 
    cplot (y, ymarker, ysize); 
    axis ('equal'); axis ('off'); title ('Original V and X'); 
    
    h_sub2 = subplot ('position', [0.05 0.1 0.25 0.3]); 
    cplot (vx, xmarker, xsize); hold on; 
    cplot (y, ymarker, ysize); 
    axis ('equal'); axis ('off'); title ('Transformed V + X');
    
    h_sub3 = subplot ('position', [0.35 0.1 0.55 0.7]);
    cplotg (vx, y, m, m_threshold); hold on;
    cMIX_plot_mixture_simple (vx, T); 
    axis ('equal'); axis ('off'); title ('Transformed V + X');
  otherwise;
end

% --- 3D data --------------------------------------------------------
else 

switch (method)
  
  case 1
    
    set (gcf, 'color', [0 0 0]);
    hold off;
    
    h_sub1 = subplot('position', [0.05 0.6 0.25 0.3]); axis ('equal'); axis ('off');
    cplot(x, xmarker, xsize); hold on; 
    cplot(y, ymarker, ysize); title ('Original V_x and Y'); 
    set(gca, 'box', 'on');
    
    h_sub2 = subplot('position', [0.05 0.1 0.25 0.3]); axis ('off');
    cplot(vx, xmarker, xsize); hold on; 
    cplot(y, ymarker, ysize); title ('Transformed V_x + Y');
    set(gca, 'box', 'on'); 
   
    h_sub3 = subplot ('position', [0.35 0.1 0.55 0.7]);
    % cMIX_plot_mixture_simple (vx, T); title ('Transformed V_x + Y');
    % keyboard
    cplotg ('black', vx, y, m, m_threshold); hold on;
    cplot (vx, xmarker, xsize); hold on; 
    cplot (y, ymarker, xsize); title ('Transformed V_x + Y'); hold on;
    axis('on'); set(gca, 'box', 'on'); rotate3d on;

    % view (6, 88);
    disp ('yahhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh');
  
  case 2
    set (gcf, 'color', [0 0 0]);
    hold off;
    
    h_sub1 = subplot ('position', [0.05 0.6 0.2 0.3]); axis ('equal');
    cplot (x, xmarker, xsize); hold on; 
    cplot (y, ymarker, ysize); title ('Original V_x and Y'); 
    set(gca, 'box', 'on');
    
    h_sub2 = subplot ('position', [0.05 0.1 0.2 0.3]);
    cplot (vx, xmarker, xsize); hold on; 
    cplot (y, ymarker, ysize); title ('Transformed V_x + Y');
    set(gca, 'box', 'on'); 

    
    h_sub3 = subplot ('position', [0.3 0.1 0.4 0.7]);
    cMIX_plot_mixture_simple (vx, T); 
    cplotg ('black', vx, y, m, m_threshold); hold on;
    cplot (vx, xmarker, xsize); hold on; 
    cplot (y, ymarker, xsize); title ('Transformed V_x + Y'); hold on;
    axis('on'); set(gca, 'box', 'on'); title ('Transformed V_x + Y');
    
    h_sub4 = subplot ('position', [0.75 0.6 0.2 0.3]); 
    % if mod(it_total,3) == 0     % this is expensive.
	
    switch (algT)
      case 'tps'
	
	cplot (y, ymarker, ysize); hold on;
	if sum(sum(c)) ~= 0
	  % ctps_plot_grid_simple('', vx1, vy, c,d); axis('equal');
	  % title ('TPS Warping');
	  vx2 = [ones(siz1,1), x] * d; vx2 = vx2(:,2:dim+1);
	  ctps_plot_grid_simple ('', x(:,1:2), y(:,1:2), c(:,1:3),d(1:3,1:3)); 

	  %cplot (vx, 'g+', xsize); hold on;
	  cplot (vx, xmarker, ysize); hold on;
	  axis('equal'); title ('TPS Warping');
	end
	
      case 'gtm_tps'
	cgtm_plot_grid_simple ('tps_style', x,y,z,w, 0); 
	axis('equal'); title ('GTM Warping');
      
      case 'gtm_gaussian'
	cgtm_plot_grid_simple ('gaussian_style', x,y,z,w, sigma_kernel); 
	axis('equal'); title ('GTM Warping');
      otherwise; disp ('no');
    end;
	
    h_sub5 = subplot ('position', [0.75 0.1 0.2 0.3]);
    vy       = m * y ./  ( (sum(m'))' * ones(1,dim));
    % [vz]     = cMIX_warp_pts (algT, z, z, c_tps, d_tps, w, sigma_kernel);
    [P,jnk] = size(z);
    
    cplot (y,  'r.', ysize); hold on;
    cplot (z,  xmarker, xsize); hold on;
    %cplot (vz, 'g+', xsize); hold on;
    %cplot (vy, ymarker, xsize); hold on;
    %cplot_2g_simple ('', z, vz, eye(P,P), 0); title ('V_z + V_y');
    
end;

end;


