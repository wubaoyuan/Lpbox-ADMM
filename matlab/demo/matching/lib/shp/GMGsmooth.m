function U = GMGsmooth(U,domainMask, numSweeps, ro, h)
% function U = GMGsmooth(U,domainMask, numSweeps, ro, h)
% ----------------------------------------------------------
% DESCRIPTION:
% 	 This function performs smoothing of the given aproximation
%
% INPUT:
% 	U                  - current aproximation
% 	domainMask  - the domain of the silhouette 
% 	h                  - grid size
% 	ro                 - right hand side of the equation
% 	numSweeps
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



% Remark: assume that U has zeros outside the domain
% so it never influences domain pixels while sweeping
% disp('GMGsmooth');
[rows, cols] = find(domainMask);
domainCoord = [rows cols];

numPixels = length(rows);
numRows = size(U,1);
numCols = size(U,2);


for currSweep = 1:numSweeps
    % here is gaus seidel relaxation
		for currPoint=1:numPixels
            
            col = domainCoord(currPoint,2);
            row = domainCoord(currPoint,1);
            
                   
            if ((col-1) == 0)
                N1=0;
            else
                N1 = U(row,col-1);
            end
            
            if ((col+1) > numCols)
                N2=0;
            else
                N2 = U(row,col+1);
            end
            
            if ((row-1) == 0)
                N3=0;
            else
                N3 = U(row-1,col);
            end
            
            if ((row+1) > numRows)
                N4=0;
            else
                N4 = U(row+1,col);
            end

            U(row,col) = -(1/4) * ...
                         (ro(row,col) * h^2 - N1 - N2 - N3 - N4);

		end
end