function [PtCs, MBs, MEs, MPs] = shpDst(siz, Pts)
% Obtain distance transform for shape data.
%
% Input
%   siz     -  size, 1 x 2
%   Pts     -  mask points, 1 x nF (cell)
%
% Output
%   PtCs    -  contour points, 1 x nF (cell), 2 x nBd
%   MBs     -  binary masks, 1 x nF (cell)
%   MEs     -  Euclidean distance transform, 1 x nF (cell)
%   MPs     -  Possion distance transform, 1 x nF (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 09-26-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
nF = length(Pts);

% per frame
[PtCs, MBs, MEs, MPs] = cellss(1, nF);
for iF = 1 : nF
    % binary mask
    MBs{iF} = maskP2M(siz, Pts{iF}([1 2], :));

    % contour
    PtC = maskCtr(MBs{iF});
    PtCs{iF} = round(PtC([2 1], :));
    MC = maskP2M(siz, PtCs{iF});

    % Euclidean distance transform
    MEs{iF} = double(bwdist(MC));

    % Possion distance transform
    MPs{iF} = GMGmain(MBs{iF});
end
