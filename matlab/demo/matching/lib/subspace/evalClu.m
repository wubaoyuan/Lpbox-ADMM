function val = evalClu(X0, H, varargin)
% Evaluate the statstics of the given cluster.
% Several statstics are available: 
%  'in':   within-class distance
%  'out':  between-class distance
%  'all':  total distance
%
% Input
%   X0      -  sample matrix, dim x n
%   H       -  indicator matrix, k x n
%   varargin
%     type  -  statstics type, {'in'} | 'out' | 'all'
%
% Output
%   val     -  value of the specific statstic
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-22-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
type = ps(varargin, 'type', 'in');

k = size(H, 1);

% within-class
ins = zeros(1, k);
for c = 1 : k
    X = cenX(X0(:, H(c, :) == 1));
    ins(c) = norm(X, 'fro');
end
in = ins * ins';

% total
X = cenX(X0);
all = norm(X, 'fro');
all = all * all;

% between-class
out = all - in;

if strcmp(type, 'in')
    val = in;
elseif strcmp(type, 'out')
    val = out;
elseif strcmp(type, 'all')
    val = all;
else
    error('unknown type');
end
