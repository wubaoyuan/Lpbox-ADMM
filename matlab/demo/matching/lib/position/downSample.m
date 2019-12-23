function [X, fPos] = downSample(X0, obj, type)
% Down samplize the original sequence.
%
% Input
%   X0      -  sequence matrix, dim x n0
%   obj     -  expected objective
%              if type = 'fps', obj is the new fps
%              if type = 'num', obj is the new frame number.
%   type    -  unit type, {'fps'} | 'num'
%
% Output
%   X       -  output matrix, dim x n
%   fPos    -  original position of each frame, 1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% default type
if ~exist('type', 'var'), type = 'fps'; end

n0 = size(X0, 2);

% objective
if strcmp(type, 'fps')
    fps0 = 120; fps = obj;
    while mod(fps0, fps) ~= 0
        fps = fps + 1;
    end

    gap = fps0 / fps;
    
elseif strcmp(type, 'num')
    maxN = obj;
    
    gap = 1;
    if n0 > maxN
        gap = ceil(n0 / maxN);
    end

else
    error(['unknown type: ' type]);
end

% down sample
fPos = 1 : gap : n0;
X = X0(:, fPos);

fprintf('down sampling of the sequence, %d frames -> %d frames\n', n0, size(X, 2));
