function    U = GMGcyc(U,domainMask,h,ro, numPreSmoothSweeps, numPostSmoothSweeps, vCycleFlag)
% function    U = GMGcyc(U,domainMask,h,ro, numPreSmoothSweeps, numPostSmoothSweeps, vCycleFlag)
% ----------------------------------------------------------
% DESCRIPTION:
% 	 This function performs one cycle of the GMG solver (recursively)
%
% INPUT:
% 	U                  - current aproximation
% 	domainMask      - the domain of the silhouette on which we should solve the poisson
% 	h                  - grid size
% 	ro                 - right hand side of the equation
% 	numPreSmoothSweeps - 
% 	numPostSmoothSweeps-
% 	vCycleFlag         - type of GMG solver (v or w)
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

% disp(strcat('GMGcyc: h=',num2str(h)));
% laplacian operator
laplacianOperator = (1/h^2) * [0 1 0; 1 -4 1; 0 1 0];
% restriction operator
restrictionOperator = (1/16) * [1 2 1; 2 4 2; 1 2 1];
% interpolation operators
interpolationOperator = (1/4) * [1 2 1; 2 4 2; 1 2 1];
numRows = size(domainMask,1);
numCols = size(domainMask,2);

ro = ro.* domainMask;% leaves zeros on the background and 
                                % the function values on the domain
% figure; imagesc(domainMask);
for cycle = 1:vCycleFlag
    
	% - presmoothing
	U = GMGsmooth(U,domainMask, numPreSmoothSweeps, ro, h);

	% - coarse grid correction
        % - compute the defect
        % the left side -  Ax vector
        % apply the laplacian operator
        Ax = conv2(U,laplacianOperator,'same');
        Ax = Ax .* domainMask;
        % |Ax-b|
        residue = ro-Ax;
        residue(~domainMask) = 0;
%         normStart = sqrt(mean(residue(domainMask).^2));
%         disp(['h: ' num2str(h) ', residue norm (start): ' num2str(normStart) ]);

        % - restrict the defect
        % restriction operator
        coarseResidue = conv2(residue,restrictionOperator,'same');
        coarseResidue = coarseResidue .* domainMask;
        
        coarseRo         = coarseResidue(1:2:end,1:2:end);
        coarseDomainMask = domainMask(1:2:end,1:2:end);
        coarseError      = zeros(size(coarseRo));
        
        % -  if the restricted defect is on a small grid
        if (nnz(coarseDomainMask) <= 5)
            % apply smoothing for direct solving
            coarseError = GMGsmooth(coarseError,coarseDomainMask, numPreSmoothSweeps, coarseRo, h*2);
        % else - apply recursion   
        else
            % call GMGcyc with the defect equation
            coarseError = GMGcyc(coarseError,coarseDomainMask,h*2,coarseRo, numPreSmoothSweeps, numPostSmoothSweeps, vCycleFlag);
        end
        % - interpolate the correction
        fineError                  = zeros(size(U));
        fineError(1:2:end,1:2:end) = coarseError;
        fineError = fineError .* domainMask;
        fineError = conv2(fineError,interpolationOperator,'same');
        fineError = fineError .* domainMask;

        % - compute the corrected aproximation
        U = U + fineError;

        % - postsmoothing
    	U = GMGsmooth(U,domainMask, numPostSmoothSweeps, ro, h);        
%                % - compute the defect
%         % the left side -  Ax vector
%         % apply the laplacian operator
%         Ax = conv2(U,laplacianOperator,'same');
%         Ax = Ax .* domainMask;
%         % |Ax-b|
%         residue = ro-Ax;
%         residue(~domainMask) = 0;
%         normStart = sqrt(mean(residue(domainMask).^2));
%         disp(['h: ' num2str(h) ', residue norm (start): ' num2str(normStart) ]);
%     figure; imagesc(residue); colorbar;
    
end