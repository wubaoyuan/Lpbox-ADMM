function flag = isY(s)
% Same to strcmpi(s, 'y').
%
% Input
%   s       -  string, 'y' | 'n'
%
% Output
%   flag    -  boolean flag
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if strcmpi(s, 'y')
    flag = true;

elseif strcmpi(s, 'n')
    flag = false;

else
    error('unknown value: %s', s);
end
