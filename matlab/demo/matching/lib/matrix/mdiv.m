function Ms = mdiv(ori, M, ks)
% Matrix division.
%
% Divide the given matrix and store the blocks the into a cell array.
%
% Example
%   input   -  ori = 'vert';
%              M = [1 2 6;
%                   3 4 5;
%                   5 6 4; 
%                   7 8 9];
%              ks = [1, 2, 1];
%   call    -  Ms = mdiv(ori, M, ks)
%   output  -  Ms{1} = [1 2 6];
%              Ms{2} = [3 4 5;
%                       5 6 4];
%              Ms{3} = [7 8 9];
%
% Input
%   ori     -  orientation, 'vert' | 'horz'
%   Ms      -  matrix, d x n
%   ks      -  #dimension, 1 x m
%
% Output
%   Ms      -  matrix set, 1 x m (cell), d_i x n_i
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-27-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m = length(ks);
Ms = cell(1, m);

if strcmp(ori, 'vert')
    head = 0;
    for i = 1 : m
        Ms{i} = M(head + 1 : head + ks(i), :);
        head = head + ks(i);
    end
    
elseif strcmp(ori, 'horz')
    head = 0;
    for i = 1 : m
        Ms{i} = M(:, head + 1 : head + ks(i));
        head = head + ks(i);
    end

else
    error('unknown direction: %s', ori);
end
