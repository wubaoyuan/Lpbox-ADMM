function Fs = vdoGray(F0s)
% Convert the color of the video from RGB to gray.
%
% Input
%   Fs      -  original video, 1 x nF (cell)
%
% Output
%   Fs      -  cropped video, 1 x nF (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
nF = length(F0s);

% convert
Fs = cell(1, nF);
for iF = 1 : nF
    Fs{iF} = imgGray(F0s{iF});
end
