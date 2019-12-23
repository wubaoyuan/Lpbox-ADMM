function anVdoF(Fs, varargin)
% Animate video data.
%
% Input
%   Fs      -  video frames, 1 x nF (cell)
%   varargin
%     save option
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% save option
[svL, path] = psSv(varargin, 'type', 'avi');

% load
if svL == 2 && exist(path, 'file')
    prIn('anVdoF', 'old');
    prOut;
    return;
end
prIn('anVdoF', 'new');

% video
nF = length(Fs);

% figure
fig = 10; figSiz = [90 160] * 2;
rows = 2; cols = 1;
Ax = iniAx(fig, rows, cols, figSiz, 'hs', [1 8]);

% text
ftSiz = 20;
s = sprintf('Frame: 0/%d', nF);
hT = shStr(s, 'ax', Ax{1}, 'ftSiz', ftSiz);
set(gca, 'Visible', 'off');

% first frame
hImg = shImg(im2double(Fs{1}), 'ax', Ax{2});

% output video
hw = vdoWIn(path, 'fps', 30);
prCIn('frame', nF, .1);
for iF = 1 : nF
    prC(iF);

    % text
    s = sprintf('Frame: %d/%d', iF, nF);
    shStrUpd(hT, s);

    % video
    shImgUpd(hImg, Fs{iF});

    hw = vdoW(hw);
end
vdoWOut(hw);
prCOut(nF);

prOut;
