function varargout = tproj(P, varargin)
% Tensor projection.
%
% Input
%   P          -  projection matrix, d x d0
%   varargin   -  original tensor, 1 x m (cell), d0 x ki x ni
%
% Output
%   varargout  -  new tensor, 1 x m (cell), d x ki x ni
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 03-20-2011
%   modify     -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m = nargin - 1;
[d, d0] = size(P);

% check
if nargout ~= m
    error('incorrect number of inputs');
end

% per tensor
for i = 1 : m
    % original tensor
    T0 = varargin{i};
    [~, k, n] = size(T0);
    
    % new tensor
    T = zeros(d, k, n);
    
    % per dimension
    for j = 1 : n
        % project
        T(:, :, j) = P * T0(:, :, j);
    end
    
    % store
    varargout{i} = T;
end
