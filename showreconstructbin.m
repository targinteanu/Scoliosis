function im = showreconstructbin(orig)

imwidth = size(orig); imheight = imwidth(1); imwidth = imwidth(2);

marker = false(size(orig));

figure; himg = imshow(orig); 
[c, r] = getpts; 

pts = (c-1)*imheight + r;
pts = floor(pts);

marker(pts) = true; 
himg.CData = marker;

im = imreconstruct(marker, orig);
%himg.CData = im;

end