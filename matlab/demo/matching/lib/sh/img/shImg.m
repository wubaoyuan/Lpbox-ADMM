function ha = shImg(F, varargin)
% Show image.
%
% Input
%   F        -  image, h x w x nChan
%   varargin
%     show option
%
% Output
%   ha
%     haImg  -  figure handle
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 01-11-2012

% show option
psSh(varargin);

% plot image
ha.haImg = imshow(F);
