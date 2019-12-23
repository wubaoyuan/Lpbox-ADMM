function ha = shImgBox(Fs, boxs, varargin)
% Plot images in boxs.
%
% Input
%   Fs       -  image, 1 x m (cell)
%   boxs     -  bounding box, 1 x m (cell), 2 x 2
%   varargin
%     show option
%
% Output
%   ha
%     hImgs  -  image handle, 1 x m
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% show option
psSh(varargin);

% dimension
m = length(Fs);

hold on;

% plot image
hImgs = cell(1, m);
for i = 1 : m
    hImgs{i} = image(boxs{i}(1, :), boxs{i}(2, :), Fs{i});
end

% colormap
colormap(gray(256));

axis ij equal;

% store
ha.hImgs = hImgs;
