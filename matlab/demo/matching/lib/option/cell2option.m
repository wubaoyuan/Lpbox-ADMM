function option = cell2option(array)
% Convert the cell (string-value) array to a struct.
% 
% Example
%   input   -  array = {'name', 'feng', 'age', 13, 'toefl', [24 25 15 26]}
%   call    -  option = cell2option(array)
%   output  -  option.name = 'feng'
%              option.age = 25
%              option.toefl = [24 25 15 26]
%
% Input
%   array   -  cell array, 1 x (2 x m) 
%
% Output
%   option  -  struct
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

m = round(length(array) / 2);

if m == 0
    option = [];
    return;
end

for i = 1 : m
    p = i * 2 - 1;
    option.(array{p}) = array{p + 1};
end
