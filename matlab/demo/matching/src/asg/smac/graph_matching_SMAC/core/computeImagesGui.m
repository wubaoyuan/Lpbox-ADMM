function [ImageOriginale , imageCouleur , image1] = computeImagesGui(image0 , p,q);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.



ImageOriginale = rescaleImage(double(image0));

[nr,nc,nb] = size(image0);
if (nb>1),
    I=double(rgb2gray(image0));
else
    I = double(image0);
end

image1 = imresize(I,[p q],'bicubic');
image1 = rescaleImage(image1);

%image1 = histeq(image1);


I1 = double(image0);
imageCouleur = imresize(I1,[p q],'bicubic');
if nb == 1
    imageCouleur = repmat(imageCouleur,[1,1,3]);
end
imageCouleur = rescaleImage(imageCouleur);

