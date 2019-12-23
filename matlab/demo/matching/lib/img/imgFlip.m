function F = imgFlip(F0, dire)
% Flip image.
%
% Input
%   F0      -  original image
%   dire    -  direction, {'horz'} | 'vert'
% 
% Output
%   F       -  new image
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-23-2012

if ~exist('dire', 'var') || strcmp(dire, 'horz')
    if ndims(F0) == 2
        F = F0(:, end : -1 : 1);
    else
        F = F0(:, end : -1 : 1, :);
    end
        
elseif strcmp(dire, 'vert')
    if ndims(F0) == 2
        F = F0(end : -1 : 1, :);
    else
        F = F0(end : -1 : 1, :, :);
    end
    
else
    error('unknown direction: %s', dire);
end
