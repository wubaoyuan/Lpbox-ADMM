function M = mcat(ori, Ms)
% Matrix concatenation.
%
% Input
%   ori     -  orientation, 'diag' | 'vert' | 'horz' | 'full'
%   Ms      -  matrix set, m1 x m2 (cell), d_i x n_i
%
% Output
%   M       -  matrix
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-27-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if strcmp(ori, 'diag')
    M = blkdiag(Ms{:});

elseif strcmp(ori, 'vert')
    M = vertcat(Ms{:});

elseif strcmp(ori, 'horz')
    M = horzcat(Ms{:});

elseif strcmp(ori, 'full')
    M = cell2mat(Ms);

else
    error('unknown direction: %s', ori);
end
