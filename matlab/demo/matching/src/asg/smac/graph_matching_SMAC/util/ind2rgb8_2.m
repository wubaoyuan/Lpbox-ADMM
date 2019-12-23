function RGB=ind2rgb8_2(X, CMAP);
% same as a function which is hidden...
% function RGB = ind2rgb8(X, CMAP)
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.


%IND2RGB8 Convert an indexed image to a uint8 RGB image.
%
%   RGB = IND2RGB8(X,CMAP) creates an RGB image of class uint8.  X must be
%   uint8, uint16, or double, and CMAP must be a valid MATLAB colormap.
%
%   Example 
%   -------
%      % Convert the 'boston.tif' image to RGB.
%      [X, cmap, R, bbox] = geotiffread('boston.tif');
%      RGB = ind2rgb8(X, cmap);
%      mapshow(RGB, R);
%
%   See also IND2RGB.

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/02/01 21:57:26 $ 


% RGB = ind2rgb8c(X, CMAP);


RGB = mex_ind2rgb8(X, CMAP);
%{
% TODO:ugly...find better solution (on super, function is opaque)
try
    RGB = ind2rgb8c(X, CMAP);
catch
    try
        dir_old=pwd;
        dir_new=fullfile(matlabroot,getpath('toolbox/images/images/private'));
        cd(dir_new);
        RGB = ind2rgb8c(X, CMAP);
        cd(dir_old);
    catch
        cd(dir_old);
        dispStack();
        error(lasterr);
    end
end
%}
