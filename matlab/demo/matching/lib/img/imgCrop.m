function [F, isIn] = imgCrop(F0, box, varargin)
% Crop image.
%
% Input
%   F0      -  original image, h0 x w0 x 3 (uint8) | h0 x w0 (double) | h0 x w0 (uint8)
%   box     -  bounding box, dim (= 2) x 2. See function boxAnd for the details about the setting of a box
%   varargin
%     sca   -  re-scaling factor (before cropping), {1}
%
% Output
%   F       -  cropped image, h x w x 3 (uint8) | h x w (double) | h x w (uint8)
%   isIn    -  position flag
%              1: Entire box is in the whole image.
%              0: Part of the box is out of the image.
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
sca = ps(varargin, 'sca', 1);

% image size
[sizF0, nC, isU] = imgInfo(F0);

% image rescaling
sizF = round(sizF0 / sca);
if sca ~= 1
    F0 = imresize(F0, sizF);
end

% box size
box0 = [1, sizF0(1); ...
        1, sizF0(2)];
siz = box(:, 2) - box(:, 1) + 1;

% box intersection
boxNew = boxAnd({box0, box});
isIn = isequal(boxNew, box);

% box position
head = boxNew(:, 1) - box(:, 1) + 1;
tail = head + boxNew(:, 2) - boxNew(:, 1);

% crop
F = imgNew(siz, nC, isU);
if nC == 1
    F(head(1) : tail(1), head(2) : tail(2)) = F0(boxNew(1, 1) : boxNew(1, 2), boxNew(2, 1) : boxNew(2, 2));
else
    F(head(1) : tail(1), head(2) : tail(2), :) = F0(boxNew(1, 1) : boxNew(1, 2), boxNew(2, 1) : boxNew(2, 2), :);
end
