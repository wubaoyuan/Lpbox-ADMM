function value = psY(option, name, default)
% Fetch the content of the specified field within a struct or string array.
% If the value is 'y', return 1, else 0.
%
% Input
%   option   -  struct or string cell array
%   name     -  field name
%   default  -  default field value
%
% Output
%   value    -  the field value
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

v = ps(option, name, default);
value = isY(v);
