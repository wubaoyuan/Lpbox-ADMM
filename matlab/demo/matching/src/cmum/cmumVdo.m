function [Fs, siz, nF] = cmumVdo(src, varargin)
% Obtain video sequence for CMU Motion source.
%
% Input
%   src     -  cmum src
%   varargin
%     save option
%
% Output
%   Fs      -  frame matrix, 1 x nF (cell)
%   siz     -  size of image, 1 x 2
%   nF      -  #frames
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-30-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-03-2013

% save option
prex = src.tag;
[svL, path] = psSv(varargin, 'fold', 'cmum/vdo', ...
                             'prex', prex, ...
                             'subx', 'vdo');

% load
if svL == 2 && exist(path, 'file')
    prInOut('cmumVdo', 'old, %s', prex);
    [Fs, siz, nF] = matFld(path, 'Fs', 'siz', 'nF');
    return;
end
prIn('cmumVdo', 'new, %s', prex);

% path
wsPath = cmumPath(src);

% video
Fs = vdoFAll(wsPath.vdo, 'comp', 'img');

% size
nF = length(Fs);
siz = imgInfo(Fs{1});

% save
if svL > 0
    save(path, 'Fs', 'siz', 'nF');
end

prOut;
