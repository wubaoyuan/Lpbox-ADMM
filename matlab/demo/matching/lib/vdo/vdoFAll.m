 function Fs = vdoFAll(fpath, varargin)
% Read all frames from a video source.
%
% Input
%   fpath   -  video path
%   varargin
%
% Output
%   Fs      -  frame set, 1 x nF (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-03-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-04-2013

% open
hr = vdoRIn(fpath, varargin{:});

% dimension
nF = hr.nF;

% read all frames
Fs = cell(1, nF);
for iF = 1 : nF
    Fs{iF} = vdoR(hr, iF);
    
    if isempty(Fs{iF})
        Fs(iF : end) = [];
        break;
    end
end

% close
vdoROut(hr);
