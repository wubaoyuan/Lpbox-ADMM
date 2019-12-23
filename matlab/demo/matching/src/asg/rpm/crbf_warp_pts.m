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


% ------------------------------------------------------------------- 
% crbf_warp_pts.m
% ------------------------------------------------------------------- 
% Warp points using RBF.
%
% Usage:
% [x1] = cgtm_warp_pts (x,z, w, sigma_kernel);
% 
% 01/26/00

function [x1] = crbf_warp_pts (x,z, w, sigma_kernel);

[phi] = crbf_gen (x, z, sigma_kernel);
x1    = phi * w;
