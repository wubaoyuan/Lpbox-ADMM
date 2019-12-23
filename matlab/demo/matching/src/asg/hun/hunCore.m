function [assignment, cost] = hunCore(D)
% Compute optimal assignment by Hungrian algorithm. 
%
% Remark
%   The distance matrix must be non-negative.
%   The distance matrix may contain infinite values (forbidden assignments).
%   Internally, the infinite values are set to a very large finite number, 
%   so that the Munkres algorithm itself works on finite-number matrices. 
%   Before returning the assignment, all assignments with infinite distance 
%   are deleted (i.e. set to zero).
%
% Input
%   D           -  distance, k1 x k2
%
% Output
%   assignment  -  index, k1 x 1
%   cost        -  cost, 1 x 1
%
% History
%   create      -  Markus Buehren (www.Lss.uni-stuttgart.de), 12-14-2004
%   modify      -  Feng Zhou (zhfe99@gmail.com), 02-28-2012

% save original distMatrix for cost computation
D0 = D;

% check for negative elements
if any(D(:) < 0)
    error('All matrix elements have to be non-negative.');
end

% dimension
[nOfRows, nOfColumns] = size(D);

% check for infinite values
finiteIndex = isfinite(D);
infiniteIndex = find(~finiteIndex);
if ~isempty(infiniteIndex)
    % set infinite values to large finite value
    maxFiniteValue = max(max(D(finiteIndex)));
    if maxFiniteValue > 0
        infValue = abs(10 * maxFiniteValue * nOfRows * nOfColumns);
    else
        infValue = 10;
    end
    
    D(infiniteIndex) = infValue;
end

% memory allocation
coveredColumns = zeros(1, nOfColumns);
coveredRows = zeros(nOfRows, 1);
starMatrix = zeros(nOfRows, nOfColumns);
primeMatrix = zeros(nOfRows, nOfColumns);

% preliminary steps
if nOfRows <= nOfColumns
    minDim = nOfRows;
	
    % find the smallest element of each row
    minVector = min(D, [], 2);
	
    % subtract the smallest element of each row from the row
    D = D - repmat(minVector, 1, nOfColumns);
	
    % Steps 1 and 2
    for row = 1:nOfRows
        for col = find(D(row,:)==0)
            if ~coveredColumns(col) %~any(starMatrix(:,col))
                starMatrix(row, col) = 1;
                coveredColumns(col)  = 1;
                break
            end
        end
    end
    
% nOfRows > nOfColumns    
else 
    minDim = nOfColumns;
    
% find the smallest element of each column
    minVector = min(D);
    
% subtract the smallest element of each column from the column
    D = D - repmat(minVector, nOfRows, 1);
    
% Steps 1 and 2
    for col = 1:nOfColumns
        for row = find(D(:,col)==0)'
            if ~coveredRows(row)
                starMatrix(row, col) = 1;
                coveredColumns(col)  = 1;
                coveredRows(row)     = 1;
                break
            end
        end
    end
    coveredRows(:) = 0; % was used auxiliary above
end

if sum(coveredColumns) == minDim
% algorithm finished
    assignment = buildassignmentvector(starMatrix);
else
% move to step 3
    [assignment, D, starMatrix, primeMatrix, ...
     coveredColumns, coveredRows] = step3(D, starMatrix, ...
                                          primeMatrix, coveredColumns, coveredRows, minDim);
end

% compute cost and remove invalid assignments
[assignment, cost] = computeassignmentcost(assignment, D0, nOfRows);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function assignment = buildassignmentvector(starMatrix)

[maxValue, assignment] = max(starMatrix, [], 2);
assignment(maxValue == 0) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, cost] = computeassignmentcost(assignment, distMatrix, nOfRows)

rowIndex   = find(assignment);
costVector = distMatrix(rowIndex + nOfRows * (assignment(rowIndex)-1));
finiteIndex = isfinite(costVector);
cost = sum(costVector(finiteIndex));
assignment(rowIndex(~finiteIndex)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step2(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim)

% cover every column containing a starred zero
maxValue = max(starMatrix);
coveredColumns(maxValue == 1) = 1;

% algorithm finished
if sum(coveredColumns) == minDim
    assignment = buildassignmentvector(starMatrix);
    
% move to step 3    
else
    [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim)

zerosFound = 1;
while zerosFound
    
    zerosFound = 0;		
    for col = find(~coveredColumns)
        for row = find(~coveredRows')
            if distMatrix(row,col) == 0
                
                primeMatrix(row, col) = 1;
                starCol = find(starMatrix(row,:));
                
                % move to step 4
                if isempty(starCol)
                    [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step4(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, row, col, minDim);
                    return;
                    
                % go on in next column
                else
                    coveredRows(row)        = 1;
                    coveredColumns(starCol) = 0;
                    zerosFound              = 1;
                    break;
                end
            end
        end
    end
end

% move to step 5
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step5(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step4(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, row, col, minDim)

newStarMatrix = starMatrix;
newStarMatrix(row,col) = 1;

starCol = col;
starRow = find(starMatrix(:, starCol));

while ~isempty(starRow)

% unstar the starred zero
    newStarMatrix(starRow, starCol) = 0;
    
% find primed zero in row
    primeRow = starRow;
    primeCol = find(primeMatrix(primeRow, :));
    
% star the primed zero
    newStarMatrix(primeRow, primeCol) = 1;
    
% find starred zero in column
    starCol = primeCol;
    starRow = find(starMatrix(:, starCol));
    
end
starMatrix = newStarMatrix;

primeMatrix(:) = 0;
coveredRows(:) = 0;

% move to step 2
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step2(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step5(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim)

% find smallest uncovered element
uncoveredRowsIndex    = find(~coveredRows');
uncoveredColumnsIndex = find(~coveredColumns);
[s, index1] = min(distMatrix(uncoveredRowsIndex,uncoveredColumnsIndex));
[s, index2] = min(s); %#ok
h = distMatrix(uncoveredRowsIndex(index1(index2)), uncoveredColumnsIndex(index2));

% add h to each covered row
index = find(coveredRows);
distMatrix(index, :) = distMatrix(index, :) + h;

% subtract h from each uncovered column
distMatrix(:, uncoveredColumnsIndex) = distMatrix(:, uncoveredColumnsIndex) - h;

% move to step 3
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);
