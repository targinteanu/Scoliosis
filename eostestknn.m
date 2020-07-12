k = 7;
fp = 'C:\Users\Toren\Desktop\scoliosis\EOS patient 7\selected DICOM\cor';
img = dicomread(fp);
X1 = repmat(1:size(img,1), size(img,2), 1);
X2 = repmat((1:size(img,2))', 1, size(img,1));
X1 = X1'; X2 = X2';
km = kmeans([double(img(:)), X1(:), X2(:)], k);
img2 = zeros(size(img)); img2(:) = km;
figure; subplot(121); imshow(img); subplot(122); imshow(label2rgb(img2));