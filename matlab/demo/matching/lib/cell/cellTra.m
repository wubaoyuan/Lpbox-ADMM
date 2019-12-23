function Ms = cellTra(M0s)
% Matrix transpose.
%
% Example
%   input   -  M0s = {[1 2 3], [4 5; 6 7]};
%   call    -  Ms = cellTra(M0s);
%   output  -  Ms = {[1; 2; 3], [4 6; 5 7]};
%
% Input
%   M0s     -  original matrix, 1 x m (cell)
%
% Output
%   Ms      -  new matrix, 1 x m (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-27-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

m = length(M0s);
Ms = M0s;
for i = 1 : m
    Ms{i} = M0s{i}';
end
