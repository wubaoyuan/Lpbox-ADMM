function posF = cutPart(nF0, maxF)
% Get a part of frames that are consecutive in order.
%
% Input
%   nF0     -  number of frames available
%   maxF    -  maximum number of frames in use
%
% Output
%   posF    -  used frame position, 1 x nF
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if nF0 < maxF
    nF = nF0 * (1 - rand(1) / 10);
else
    nF = maxF * (1 - rand(1) / 10);
end
nF = floor(nF);

head = float2block(rand(1), nF0 - nF + 1);
posF = head : head + nF - 1;
