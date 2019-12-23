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


% cICP_findneighbours.m
% ------------------------------------------------------------------- 
% Find nearest neighbour in template t for each point in x. 3d/2d
%
% Usage: [M, distances] = cICP_findneighbours (x, t);
% Notes: 
%         M -- M(1) = 10, y(10) is nearest from x(1);
%         distance -- distance (x(i)-t(j))^2
% ------------------------------------------------------------------- 
% Last modified: 02/07/98

function [M, distances] = cICP_findneighbours (x, t);

[m, dim] = size(x);
[n, dim] = size(t);
M = zeros (m,1);


% |x-t| matrices:
% ---------------
xttmp = zeros (n, m);
for i=1:dim
  xtmp = ones(n,1) * x(:,i)';
  ttmp = t(:,i)  * ones(1,m); 
  xttmp = xttmp + (xtmp - ttmp) .* (xtmp - ttmp);
end;



% M + min_dist list:
% ------------------
[min_dist, min_index] = min(xttmp);
distances = (sqrt(min_dist))';
M         = min_index';












% theta = transformation (1);
% t(1)  = transformation (2);
% t(2)  = transformation (3);
% s     = transformation (4);

% cs = cos(theta); sn = sin(theta); R = [cs, -sn; sn, cs];
% ry = y*R'*s;
% tr_y = ry;
% for k=1:2, tr_y(:,k) = tr_y(:,k) + t(k);, end

% for i=1:xmax
%
%   dist_tmp = 0;
%   dist_min = 0;
%   index    = -1;
%   
%   dist_min = ((tr_y(1,:) - x(i,:)) * (tr_y(1,:) - x(i,:))');
%   index    = 1;
%   
%   for j=2:ymax
% 
%     dist_tmp = ((tr_y(j,:) - x(i,:)) * (tr_y(j,:) - x(i,:))');
%    if dist_tmp < dist_min
%       dist_min = dist_tmp;
%       index = j;
%     end;
%   end;
%   
%   M(i) = index;
%   distances(i) = sqrt(dist_min);
% end;
