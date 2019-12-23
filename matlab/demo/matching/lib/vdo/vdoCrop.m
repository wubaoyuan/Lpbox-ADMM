function [Fs, isIns] = vdoCrop(F0s, Box, varargin)
% Crop video.
%
% Input
%   F0s     -  original video, 1 x nF (cell)
%   Box     -  bounding box, 2 x 2 x nF. See function boxAnd for the details about the setting of a box
%   varargin
%     sca   -  re-scaling factor, {1}
%
% Output
%   Fs      -  cropped video, 1 x nF (cell)
%   isIns   -  position flag, 1 x nF
%              1: Entire box is in the whole video.
%              0: Part of the box is out of the video.
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
sca = ps(varargin, 'sca', 1);

% dimension
nF = length(F0s);

% check box size
sizs = squeeze(Box(:, 2, :) - Box(:, 1, :)) + 1;
siz = sizs(:, 1);
if any(sizs(1, :) ~= siz(1)) || any(sizs(2, :) ~= siz(2))
    error('incorrect bounding box');
end

% crop
Fs = cell(1, nF);
isIns = zeros(1, nF);
for iF = 1 : nF
    [Fs{iF}, isIns(iF)] = imgCrop(F0s{iF}, Box(:, :, iF), 'sca', sca);
end
