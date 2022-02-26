img = imread('2021-08-18 (4).png');

img = double(img)/255;
imgred = img(:,:,1)==1;
N = [imgred(2:end,:); false(size(imgred(1,:)))];
S = [false(size(imgred(1,:))); imgred(1:(end-1),:)];
W = [imgred(:,2:end), false(size(imgred(:,1)))];
E = [false(size(imgred(:,1))), imgred(:,1:(end-1))];
NW = [imgred(2:end,2:end), false(size(imgred(2:end,1))); false(size(imgred(1,:)))];
SE = [false(size(imgred(1,:))); false(size(imgred(2:end,1))), imgred(1:(end-1),1:(end-1))];
NE = [false(size(imgred(2:end,1))), imgred(2:end,1:(end-1)); false(size(imgred(1,:)))];
SW = [false(size(imgred(1,:))); imgred(1:(end-1),2:end), false(size(imgred(2:end,1)))];


outln_ = (N|S|W|E|NW|SE|NE|NW) & (~imgred);
outln = false(size(img));
outln(:,:,1) = outln_; 
outln(:,:,2) = outln_; 
outln(:,:,3) = outln_;

img(outln) = 1;
figure; imshow(img)