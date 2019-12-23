function Urelaxed = relaxAfterSolvingU(U, numSweeps)
% function Urelaxed = relaxAfterSolvingU(U, numSweeps)
% ----------------------------------------------------------
% DESCRIPTION:
%   relaxes the solution  to the poisson equation with 0's outside the
%   sihlouette and -1's inside the sihlouette. This action will "blur" the
%   values of U from inside to the outside, while hardly changing the
%   inside values
%
% INPUT:
%   U - approximation of the solution
%   numSweeps  - the number of relaxation sweeps to perform (default = 4)
%
% ASSUMPTIONS:
%
% OUTPUT:
%   Urelaxed - the relaxed version of U
%   
% AFFECTED MAT FILES:
%
% ----------------------------------------------------------

chi = (U > 0);
numRows = size(U,1);
numCols = size(U,2);
if ~exist('numSweeps','var')
    numSweeps = 4;
end

% figure;
for currSweep = 1:numSweeps
    % here is gaus seidel relaxation
    for row = 1:numRows
        for col = 1:numCols

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

            % different relaxation inside and outside
            if (chi(row, col))
                U(row,col) = (1 + N1 + N2 + N3 + N4)/4;
            else
                U(row,col) = (N1 + N2 + N3 + N4)/4;
            end
        end
    end
end
Urelaxed = U;
return