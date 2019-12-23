function hw = vdoWIn(fpath, varargin)
% Create a video handler for writing.
%
% Input
%   fpath    -  video path
%   varargin
%     comp   -  compressing mode, 'none' | {'vdo'} | 'img' | 'mcv'
%     nF     -  only used when comp == 'img'
%     siz    -  only used when comp == 'img' or 'mcv', [height, width]
%     pathform  -  only used when comp == 'img'
%     idx    -  only used when comp == 'img'
%     fps    -  frame per second, {30} | ...
%
% Output
%   hw
%     vdo    -  video handler
%     comp   -  compressing mode
%     fpath  -  video path
%     siz    -  size, [height, width]
%     fps    -  frame per second
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
comp = ps(varargin, 'comp', 'vdo');
fps = ps(varargin, 'fps', 30);

% non-compressed
if strcmp(comp, 'none')
    hw.vdo = avifile(fpath, 'fps', fps);
    
    prIn('vdoWIn', '%s', fpath);
    pr('fps %d, comp %s', fps, comp);

% VideoReader
elseif strcmp(comp, 'vdo')
    hw.vdo = VideoWriter(fpath, 'Motion JPEG AVI');
    hw.vdo.FrameRate = fps;
    open(hw.vdo);
    
    prIn('vdoWIn', '%s', fpath);
    pr('fps %d, comp %s', fps, comp);

% image
elseif strcmp(comp, 'img')
    if length(fpath) > 3 && (strcmp(fpath(end - 3 : end), '.avi') || strcmp(fpath(end - 3 : end), '.m4v'))
        fpath(end - 3 : end) = [];
        hw.fpath = fpath;
    end
    if ~exist(fpath, 'dir')
        mkdir(fpath);
    end

    % info
    infopath = sprintf('%s/info.mat', fpath);
    nF = ps(varargin, 'nF', 1);
    siz = ps(varargin, 'siz', []);
    pathform = ps(varargin, 'pathform', '%.4d.jpg');
    idx = ps(varargin, 'idx', 1 : nF);
    save(infopath, 'nF', 'siz', 'fps', 'pathform', 'idx');

    hw.iF = 0;
    hw.nF = nF;
    hw.siz = siz;
    hw.pathform = pathform;
    hw.idx = idx;

    prIn('vdoWIn', '%s', fpath);
    pr('nF %d, siz %d x %d, fps %d, comp %s, pathform %s', nF, siz(1), ...
       siz(2), fps, comp, pathform);
    
% mat
elseif strcmp(comp, 'mat')

% OpenCV
elseif strcmp(comp, 'mcv')
    siz = ps(varargin, 'siz', []);
    hw.vdo = mcvVWriter(fpath, siz(2), siz(1), fps);
    
    prIn('vdoWIn', '%s', fpath);
    pr('siz %d x %d, fps %d, comp %s', siz(1), siz(2), fps, comp);
    
else
    error('unknown comp: %s', comp);
end

% store
hw.comp = comp;
hw.fpath = fpath;
hw.fps = fps;

prOut;
