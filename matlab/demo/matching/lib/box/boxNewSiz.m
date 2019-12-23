function box = boxNewSiz(box0, sca)
% Readjust the box size.
% 
% Input
%   box0    -  original box, 2 x 2
%   sca     -  scale factor, 2 x 1
%
% Output
%   box     -  new box, 2 x 2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-13-2012

mi0 = box0(:, 1);
ma0 = box0(:, 2);
siz0 = ma0 - mi0;
me0 = (mi0 + ma0) / 2;

siz = siz0 .* sca;
box = [me0 - siz ./ 2, me0 + siz ./ 2];
