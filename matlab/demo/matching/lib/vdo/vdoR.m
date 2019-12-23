function F = vdoR(hr, iF)
% Retrieve frame from video.
%
% Input
%   hr      -  video handler
%   iF      -  frame index
%
% Output
%   F       -  frame, h x w x 3 (uint8) | h x w (double) | h x w (uint8)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% VideoReader
if strcmp(hr.comp, 'vdo')
    F = read(hr.vdo, iF);

elseif strcmp(hr.comp, 'img')
    imgpath = sprintf(['%s/' hr.pathform], hr.fpath, hr.idx(iF));
    F = imread(imgpath);
    
elseif strcmp(hr.comp, 'mat')
    F = hr.Fs{iF};
    
% OpenCV
elseif strcmp(hr.comp, 'mcv')
    F = mcvVReader(hr.vdo);
end

% uint8 -> double
if strcmp(hr.form, 'double')
    F = im2double(F);
end

% rgb -> gray
if strcmp(hr.cl, 'gray')
    F = rgb2gray(F);
end
