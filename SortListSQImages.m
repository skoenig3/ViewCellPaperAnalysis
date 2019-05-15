%% Automatically Process images from Flickr Pics
% ***Should only run if get_Flickr_Seth is done running.***
% [1] Remove images that are too small or are grayscale
% [2] Resize images
% [3] Select for "interesting" Pictures
%% [1] Remove images that are too small or are grayscale
%moves small and grayscale images to seperate foriginaler. Does not delete them.
%Must confirm everything manually.

imageX = 800;% desired horizontal image size
imageY = 600;% desired vertical image size

imgdir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\FlickrImages\';
cd(imgdir)
smalldir = 'Images that are too small\';
mkdir([imgdir smalldir]);
graydir = 'Gray scale images\';
mkdir([imgdir graydir]);
a = ls;

toosmallimages = [];
grayimages = [];
for i = 1:size(a,1);
    index = strfind(a(i,:),'bmp');
    if isempty(index);
        index = strfind(a(i,:),'jpg');
        if isempty(index)
            index = strfind(a(i,:),'jpeg');
        end
    end
    if ~isempty(index);
        img = imread(a(i,:));
        if size(img,3) == 1
            grayimages= [grayimages i];
        elseif all(all(img(:,:,1) == img(:,:,2))) ...
                || all(all(img(:,:,2) == img(:,:,3))) || all(all(img(:,:,1) == img(:,:,3)))
            grayimages= [grayimages i];
        elseif size(img,1) < imageY || size(img,2) < imageX
            toosmallimages = [toosmallimages i];
        end
    end
end
%%
for ti = 1:length(toosmallimages);
    movefile(a(toosmallimages(ti),:),[imgdir smalldir])
end
for gi = 1:length(grayimages);
    movefile(a(grayimages(gi),:),[imgdir graydir]);
end

%% [2] Resize images
% assumes no .mat files in the foriginaler
% resizes images to desired size using imresize
clear, clc

imageX = 800;% desired horizontal image size
imageY = 600;% desired vertical image size

imgdir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\FlickrImages\';
resized_dir = [imgdir 'Resized\'];
mkdir(resized_dir);

a = ls;
for aa = 1:size(a,1);
    if ~isempty(strfind(a(aa,:),'bmp')) || ~isempty(strfind(a(aa,:),'jpg')) || ~isempty(strfind(a(aa,:),'_b'))...
            ||  ~isempty(strfind(a(aa,:),'jpeg'))
        img = imread(a(aa,:));
        if all(size(img) == [imageY,imageX,3]);
            imwrite(img,[resized_dir a(aa,1:end-4) '.bmp'],'bmp');
            delete([imgdir a(aa,:)]);
        else
            img = imresize(img,[imageY,imageX]);
            imwrite(img,[resized_dir a(aa,1:end-4) '.bmp'],'bmp');
            delete([imgdir a(aa,:)]);
        end
    end
end
%% [3] select for intereseting pictures
% use variability in pixel intensities (image entropy) and percent edginess 
% (high freqency content) to weed out images lacking in dynamic content.
% This is an objective process but isn't totally perfect. Essentially code
% removes pictures with a lot of background and low range of image
% intensities (i.e. images that are too bright or too dark)

clear, clc
imgdir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\FlickrImages\Resized\';
cd(imgdir)
a = ls;

% sobel filters detect edges!
sobelx = [1     2   1;
    0     0   0;
    -1    -2  -1;];

sobely = [1     2   1;
    0     0   0;
    -1    -2  -1;];

entropyvalues = NaN(1,size(a,1));
edgevalues = NaN(1,size(a,1));
for i = 1:size(a,1);
    index = strfind(a(i,:),'bmp');
    if isempty(index);
        index = strfind(a(i,:),'jpg');
        if isempty(index)
            index = strfind(a(i,:),'jpeg');
        end
    end
    if ~isempty(index);
        img = imread(a(i,:));
        img = rgb2gray(img);
        entropyvalues(i) = entropy(img);%pixel intesnity entropy
        xedges = imfilter(img,sobelx);
        yedges = imfilter(img,sobely);
        edgevalues(i) = mean2(xedges+yedges); %edgineess
    end
end

gonogo = zeros(1,size(a,1));
for i = 1:size(a,1);
    if entropyvalues(i) > 7 && (edgevalues(i) > 15) &&  (edgevalues(i) <85)
        gonogo(i) = 1;
    end
end

gooddir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\FlickrImages\Good Images\';
mkdir(gooddir);


for i = 1:size(a,1);
    if gonogo(i) == 1
        movefile(a(i,:),gooddir);
    end
end