function vdoWOut(hw)
% Close a video handler for writing.
%
% Input
%   hw      -  handler
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% non-compressed
if strcmp(hw.comp, 'none')
    hw.vdo = close(hw.vdo);

% VideoReader
elseif strcmp(hw.comp, 'vdo')
    close(hw.vdo);

elseif strcmp(hw.comp, 'img')
    
elseif strcmp(hw.comp, 'img')
    
% OpenCV    
elseif strcmp(hw.comp, 'mcv')
    mcvVWriter(hw.vdo);
    
else
    error('unknown comp: %s', hw.comp);
end

