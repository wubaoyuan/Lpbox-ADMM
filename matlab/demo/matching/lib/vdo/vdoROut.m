function vdoROut(hr)
% Close a video handler for reading.
%
% Input
%   hr      -  handler
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if strcmp(hr.comp, 'vdo')
    
elseif strcmp(hr.comp, 'img')
    
elseif strcmp(hr.comp, 'mat')

elseif strcmp(hr.comp, 'mcv')
    mcvVReader(hr.vdo);
end

