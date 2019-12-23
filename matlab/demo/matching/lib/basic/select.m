function value = select(flag, v1, v2)
% Implementation of : operator in c++.
%
% Input
%   flag    -  binary flag
%   v1      -  the first value
%   v2      -  the second value
%
% Output
%   value   -  value = flag ? v1 : v2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if flag
    value = v1;
else
    value = v2;
end
