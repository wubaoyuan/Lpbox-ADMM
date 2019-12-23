function tokens = tokeniseSlow(s, delim)
% Split the given string into separate tokens.
%
% Input
%   s       -  string, 1 x len (char)
%   delim   -  delims, 1 x m (char)
%
% Output
%   tokens  -  token list, 1 x n (cell)
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 12-30-2008
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

maxN = 10000;
len = length(s);
tokens = cell(1, maxN);

% scan the string char by char
n = 0; p = 1;
while p <= len
    
    % skip the heading empty space
    while p <= len && isContain(delim, s(p))
        p = p + 1;
    end
    
    if p > len, break; end
    
    head = p;
    while p <= len && ~isContain(delim, s(p))
        p = p + 1;
    end
    n = n + 1;
    tokens{n} = s(head : p - 1);
end

tokens(n + 1 : end) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flag = isContain(s, a)
% Return whether the character is contained in the string.
% Input
%   s     -  string, 1 x m
%   a     -  char
% Output
%   flag  -  boolean flag

for i = 1 : length(s)
    if s(i) == a
        flag = true;
        return;
    end
end

flag = false;
