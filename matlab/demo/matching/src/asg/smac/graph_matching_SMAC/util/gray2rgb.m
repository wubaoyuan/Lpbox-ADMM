function image2 = gray2rgb(image,colorMap,isRescale);
% Timothee Cour, 21-Apr-2008 17:31:23
% This software is made publicly for research use only.
% It may be modified and redistributed under the terms of the GNU General Public License.

if isempty(image)
    image2=[];
    return;
end
if nargin<3
    isRescale=1;
end

image=double(image);
if nargin<2
    colorMap=jet(128);%TODO:unused...
end
%image2 = gray2rgb(image,jet); for 'default' colormap
[p,q,r] = size(image);
image = full(image);
if r==1
%     image2=repmat(image,[1,1,3]);
    if nargin < 2
        if isRescale
            image=rescaleImage(image);
        end
        image2=repmat(image,[1,1,3]);
    else
        k=size(colorMap,1);
        
        if isRescale
%             image=rescaleImage(image);
%             image=round(rescaleImage(image,[1,k]));
            image=rescaleImage(image,[1,k]);
        else
%             image=round((k-1)*image+1);
            image=(k-1)*image+1;
        end       
%         image=uint8(image);
%         image2 = ind2rgb(image,colorMap);
%         image2=uint8(image2);
%         image2 = ind2rgb8(image, colorMap);
        image2 = ind2rgb8_2(image, colorMap);
    end
else
    image2 = image;
end

% image2=rescaleImage(image2);