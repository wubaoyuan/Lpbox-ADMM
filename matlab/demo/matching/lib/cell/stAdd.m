function a = stAdd(a, varargin)
% Add field to a given struct.
%
% Example
%   input     -  a.name = 'feng';
%   call      -  a = stFld(a, 'age', 25, 'gender', 'male');
%   input     -  a.name = 'feng';
%                a.age = 25;
%                a.gender = 'male'.
%
% Input
%   a         -  struct
%   varargin  -  field name list, 1 x m (cell)
%
% Output
%   a         -  struct      
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 02-01-2009
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

m = round(length(varargin) / 2);
for i = 1 : m
    a.(varargin{(i - 1) * 2 + 1}) = varargin{i * 2};
end
