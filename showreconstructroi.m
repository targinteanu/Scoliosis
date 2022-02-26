function im = showreconstructroi(orig)

imwidth = size(orig); imheight = imwidth(1); imwidth = imwidth(2);

figure; %himg = imshow(orig); 
marker = double(roipoly(orig));

im = imreconstruct(marker, orig);
%himg.CData = im;
imshow(im);

end