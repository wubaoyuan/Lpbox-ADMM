function boxs = imgBox(Fs, varargin)
% Obtain bounding box of image.
% 
% Input
%   Fs      -  image, 1 x m (cell)
%   varargin
%
% Output
%   boxs    -  bounding box, 1 x m (cell), 2 x 2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-02-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m = length(Fs);

boxs = cell(1, m);
for i = 1 : m
    siz = imgInfo(Fs{i});
    boxs{i} = [0, 0; siz([2 1])]';
end
