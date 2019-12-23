function pr(form, varargin)
% Prompt the information specified in the parameters.
%
% Input
%   form      -  format
%   varargin  -  object list
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 01-29-2009
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% set by function "prSet.m"
global lPr lMaPr;

if lPr <= lMaPr
    for l = 1 : lPr
        fprintf('-');
    end

    fprintf(form, varargin{:});
    fprintf('\n');
end
