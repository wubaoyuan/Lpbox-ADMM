function hw = vdoW(hw0)
% Write one frame into video.
%
% Input
%   hw0     -  original video handler
%
% Output
%   hw      -  new video handler
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 11-24-2012

hw = hw0;

% non-compressed
if strcmp(hw.comp, 'none')
    drawnow; 
    hw.vdo = addframe(hw.vdo, gcf);

% VideoReader
elseif strcmp(hw.comp, 'vdo')
    drawnow; 
    F = vdoGrab(gcf);
    writeVideo(hw.vdo, F);

elseif strcmp(hw.comp, 'img')
    hw.iF = hw.iF + 1;

    imgpath = sprintf(['%s/' hw.pathform], hw.fpath, hw.idx(hw.iF));
    drawnow;
    F = vdoGrab(gcf);
    imwrite(F.cdata, imgpath);
    
elseif strcmp(hw.comp, 'mat')
    
% OpenCV    
elseif strcmp(hw.comp, 'mcv')
    drawnow; 
    F = vdoGrab(gcf);
    mcvVWriter(hw.vdo, F.cdata);
end
