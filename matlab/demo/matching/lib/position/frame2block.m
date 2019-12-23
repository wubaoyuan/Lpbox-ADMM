function pos = frame2block(pos0, st)
% Convert frame-based position to block-based position.
%
% Input
%   pos0    -  frame-based position, 1 x n
%   st      -  stating position of segments, 1 x (m + 1)
%
% Output
%   pos     -  block-based position, 1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

n = length(pos0);

if isempty(st)
    pos = ones(1, n);
    return;
end

pos = pos0;

% position of block
j = 1;

% position of pos0
for i = 1 : n
    while j < length(st) && pos0(i) >= st(j + 1)
        j = j + 1;
    end

    if j == length(st)
        error('out of index');
    end

    pos(i) = j;
end
