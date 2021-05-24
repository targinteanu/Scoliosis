fn = 'C:\Users\Toren\Desktop\scoliosis\EOS patient 36\selected DICOM\cor';
imcor = dicomread(fn);
imcor = double(imcor); imcor = imcor - min(imcor(:)); imcor = imcor/max(imcor(:));
ifo = dicominfo(fn);
xscl = ifo.PixelSpacing(1); zscl = ifo.PixelSpacing(2);

%[outln, vertx, verty] = roipoly(imcor); 
load('EOStestROI.mat');

% img is a double img [0, 1]; outln is a BW outline

% filtering 
im = imreconstruct(double(outln), imcor); 
%imgFiltered = imcor - im; 
imgFiltered = im;
imgFiltered = imgFiltered .* outln;
imgFiltered = imgFiltered/max(imgFiltered(:));

figure; 
subplot(131); imshow(imcor, []); hold on; visboundaries(outln); 
subplot(132); imshow(im, []);
subplot(133); imshow(imgFiltered, []);

rcref = round(CenterOfMass(outln));

%% various processes 
imcorfill = imfill(imcor);

ws = label2rgb(watershed(imcorfill)); 
ff = labeloverlay(imcorfill, grayconnected(imcorfill, rcref(1), rcref(2)));
fm = labeloverlay(imcorfill, imsegfmm(graydiffweight(imcorfill, outln, 'GrayDifferenceCutoff', 25), outln, graythresh(imcorfill)));
km = label2rgb(imsegkmeans(single(imcorfill), 10));
ac = label2rgb(activecontour(imcorfill, outln));

ej_sob = edge(imcorfill, 'Sobel');
ej_prw = edge(imcorfill, 'Prewitt');
ej_rob = edge(imcorfill, 'Roberts');
ej_log = edge(imcorfill, 'log');
ej_can = edge(imcorfill, 'Canny');

[circC, circR, circM] = imfindcircles(imcorfill, [130,170]);
[Gmag, Gdir] = imgradient(imcorfill);

figure; 
subplot(331); imshow(ws); title('watershed'); 
subplot(332); imshow(ff); title('flood fill'); 
subplot(333); imshow(fm); title('fast marching'); 
subplot(334); imshow(km); title('k means');
subplot(335); imshow(ac); title('active contour');
subplot(336); imshow(ej_sob); title('edge'); 
subplot(337); imshow(Gmag); title('grad mag'); 
subplot(338); imshow(Gdir); title('grad dir'); 
subplot(339); imshow(imcorfill); hold on; viscircles(circC, circR);

figure; 
subplot(321); imshow(ej_sob); title('Sobel');
subplot(322); imshow(ej_prw); title('Prewitt');
subplot(323); imshow(ej_rob); title('Roberts');
subplot(324); imshow(ej_log); title('log');
subplot(325); imshow(ej_can); title('Canny');