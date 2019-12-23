function [i1,i2,j1,j2] = position2lim(left, bottom, width, height)
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if nargin == 1
    bottom=left(2);
    width = left(3);
    height = left(4);
    left = left(1);
end
pos0 = get(0,'ScreenSize');
height0 = pos0(4);
% [left0, bottom0, width0, height0] = get(0,'ScreenSize');

i2 = height0 - bottom+1;
j1 = left;
i1 = i2 - height+1;%bottom - height + 1;
j2 = left + width - 1;

if nargout == 1
    i1 = [i1,i2,j1,j2];
end