function Ms = cellVec(M0s)
% Matrix vectorization.
%
% Example
%   input   -  M0s = {[1 2 3], [4 5; 6 7]};
%   call    -  Ms = cellVec(M0s);
%   output  -  Ms = {[1; 2; 3], [4; 6; 5; 7]};
%
% Input
%   M0s     -  original matrix, 1 x m (cell)
%
% Output
%   Ms      -  new matrix, 1 x m (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 09-25-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

m = length(M0s);
Ms = cell(1, m);
for i = 1 : m
    Ms{i} = M0s{i}(:);
end
