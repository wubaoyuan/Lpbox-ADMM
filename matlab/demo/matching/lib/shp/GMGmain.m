function U = GMGmain(domainMask)
% function U = GMGmain(domainMask)
% ----------------------------------------------------------
% DESCRIPTION:
% 	 This function gets as an input a logical matrix that represents the
% 	 sihlouette (zeros outside and ones inside) and returns matrix that 
%    represents an aproximate solution U  to the poisson equation 
%    Uxx+Uyy=-1 defined on the domain of the silhouette with zero boundary
% 	 conditions at the contour of the sihlouette
%
% INPUT:
%   domainMask - logic matrix that holds the sihlouete domain
%
% 
% OUTPUT:
%   U - approximate solution to the poisson equation
%
% AUTHOR:
%   Lena Gorelick
%
% YEAR:
%   2003
%
% INSTITUTION:
%   Weizmann
%
% ----------------------------------------------------------

% disp ('GMGmain');

    
    % first aproximation
    U=zeros(size(domainMask),'single');
    % right hand side of the equation
    % IS MEANINGFUL ONLY IN THE domainMask
    ro=-1*single(domainMask);
   
    %-------------------------------------
    % S Y S T E M      P A R A M E T E R S
    %-------------------------------------
    % the grid size
    h=1;
    numPreSmoothSweeps  = 2;
    numPostSmoothSweeps = 2;
    % type of GMG
    vCycleFlag = 1;
    numCycles = 1;
    for i=1:numCycles
        U = GMGcyc(U,domainMask,h,ro, numPreSmoothSweeps, numPostSmoothSweeps, vCycleFlag);
    end
return
