function prIn(nm, form, varargin)
% % Start a propmter for displaying information.
%
% Input
%   nm        -  name
%   form      -  format
%   varargin  -  object list
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 01-29-2009
%   modify    -  Feng Zhou (zhfe99@gmail.com), 03-02-2012

% variables set in "prSet.m"
global lPr;

% check
if nargin == 0 || ~ischar(nm)
    error('incorrect input for prIn');
end

% init
if isempty(lPr)
    prSet(3);
end

% print
if nargin == 1
    pr('%s', nm);
else
    pr(['%s: ' form], nm, varargin{:});
end

% insert
lPr = lPr + 1;
