% Robust Point Matching (RPM) Demo (version 20000427):
% ----------------------------------------------------
% Copyright (C) 2000 Haili Chui, Anand Rangarajan
% 
% Authors: Haili Chui and Anand Rangarajan
% Date:    04/27/2000
% 
% Contact Information:
%
% Haili Chui:		chui@noodle.med.yale.edu
% Anand Rangarajan:	anand@noodle.med.yale.edu
% 
% Terms:	  
% 
% The source code (M-files) are provided under the
% terms of the GNU General Public License with an explicit
% clause permitting the execution of the M-files from within
% a MATLAB environment. See the LICENSE file for details.
%
%


% ---------------------------------------------------------------
% ctps_plot_grid_gen.m
% ---------------------------------------------------------------
% Generate grid points for displaying TPS deformation.
%
% Usage:
% [grid_pts] = ctps_plot_grid_gen (x);
% [grid_pts] = ctps_plot_grid_gen (x, resolution, resolution_grid);
% 
% 01/26/00

function [grid_pts] = ctps_plot_grid_gen (x, resolution, resolution_grid);

% check input:
% ------------
if (nargin == 1);      % input (x), set the other 2.
  resolution = 4;
  resolution_grid = 3;
elseif (nargin == 3);  % input (x, resolution, resolution_grid);
  ;
else
  disp ('# ERROR #: ctps_plot_grid_gen -- wrong input!');
  help ctps_plot_grid_gen;
end;



% set grid range:
% ---------------
x1=min(x(:,1)); x2=max(x(:,1));
y1=min(x(:,1)); y2=max(x(:,1));

% expand a little bit:
expand_ratio = 5;
x1 = x1 - (x2-x1)/expand_ratio;
x2 = x2 + (x2-x1)/expand_ratio;
y1 = y1 - (y2-y1)/expand_ratio;
y2 = y2 + (y2-y1)/expand_ratio;

% square grid:
a = min(x2-x1, y2-y1) / resolution;
grid_step = a / resolution_grid;

rows = ceil((yrange(2)-yrange(1)) / a + 1);
cols = ceil((xrange(2)-xrange(1)) / a + 1);

y2 = y1 + (rows-1)*resolution_grid*grid_step;
x2 = x1 + (cols-1)*resolution_grid*grid_step;

% generate grid points:
grid_pts = [];


% row by row first: 
xstart = x1; xend   = x2; 
x_step = grid_stp;
y_step = (y2-y1);

xnow = xstart;
ynow = ystart;

% 4 direction: up, left, down, left;
dlist = [0 y_step; x_step 0; 0 -y_step; x_step, 0];
while (xnow < xend)
  for i=1:4
    xnext = xnow + dlist(i,1);
    ynext = ynow + dlist(i,1);
    
    tmp = [[xnow:grid_step:xnext]',... 
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % stopped here !!!!
      % not working yet.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% points_row = floor( (xrange(2)-xrange(1))/grid_step ); 
points_row = (cols-1) * resolution_grid + 1; % two ending points.
for i=1:rows
  tmp_row = [[xrange(1):grid_step:xrange(2)]', ...
	ones(points_row,1) * (i-1) * a + yrange(1)];
  grid_pts = [grid_pts; tmp_row];
end;


% points_col = floor((yrange(2)-yrange(1))/grid_step ); 
points_col = (rows-1) * resolution_grid + 1; % two ending points.
for j=1:cols
  tmp_col = [ones(points_col,1) * (j-1) * a + xrange(1), ...
	[yrange(1):grid_step:yrange(2)]'];
  grid_pts = [grid_pts; tmp_col];
end;

