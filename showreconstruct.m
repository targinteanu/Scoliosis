function im = showreconstructroi(orig)

imwidth = size(orig); imheight = imwidth(1); imwidth = imwidth(2);

marker = zeros(size(orig));

figure; himg = imshow(orig); 
[c, r] = getpts; 

pts = (c-1)*imheight + r;
pts = floor(pts);

marker(pts) = 1; 

im = imreconstruct(marker, orig);
himg.CData = im;

end