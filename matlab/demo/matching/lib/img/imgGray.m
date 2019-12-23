function F = imgGray(F0)
% Convert the color of the image from RGB to gray.
%
% Input
%   F0      -  original image
% 
% Output
%   F       -  new image
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 06-02-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if isa(F0, 'uint8')
    F = im2double(F0);
end

F = rgb2gray(F);
