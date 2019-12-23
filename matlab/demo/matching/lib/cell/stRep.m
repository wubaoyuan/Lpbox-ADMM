function A = stRep(a, varargin)
% Replicate a struct into a cell matrix.
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
%   A         -  struct matrix
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 04-07-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

A = cell(varargin{:});
for i = 1 : numel(A)
    A{i} = a;
end
