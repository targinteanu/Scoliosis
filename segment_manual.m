%% load volume 
%filepath = uigetdir; 
filepath = 'C:\Users\Toren\Desktop\scoliosis\Patient 22\Lumbar Spine C- Ct\L-SPINE - 2';

vol = dicomreadvol(filepath); 

slc = 175;

%% basic filtering 
figure; imshow(vol(:,:,slc))
marker = zeros(size(vol)); %marker(305, 141, slice) = 1; 
marker(270,141,slc) = 1;
bg = imreconstruct(marker, vol); 
volfilt = vol - bg; 
figure; imshow(volfilt(:,:,slc));
volhist = histeq(volfilt);
volhist = volhist - min(volhist(:)); volhist = volhist/max(volhist(:));
figure; imshow(volhist(:,:,slc));
%% thresholding 
volthreshed = volhist > .6; 
figure; imshow(volthreshed(:,:,slc))
%% erosion/dilation 
%{
se = strel('sphere', 5);
volerode = imerode(volthreshed, se);
%figure; imshow(volerode(:,:,slice));
voldil = imdilate(volerode, se);
figure; imshow(voldil(:,:,slice));
%}
%% CC & filtering based on size 
%CC = bwconncomp(volerode, 6); 
CC = bwconncomp(volthreshed, 6);
L = labelmatrix(CC);
figure; imshow(label2rgb(L(:,:,slc)));

objs = CC.PixelIdxList;
minsize = 10000; 
szs = zeros(size(objs));
for n = 1:length(szs)
    szs(n) = length(objs{n});
end
objs = objs(szs > minsize); % objs{n} is an index matrix, not volume 
szs = szs(szs > minsize);

shapes = cell(size(objs)); Z = false(size(volthreshed));
for n = 1:length(objs)
    I = Z; I(objs{n}) = true;
    shapes{n} = I; % shapes{n} is a volume matrix
end

CCfilt = CC;
CCfilt.PixelIdxList = objs;
CCfilt.NumObjects = length(objs);
Lfilt = labelmatrix(CCfilt);
figure; imshow(label2rgb(Lfilt(:,:,slc)));

%% watershed 
[~,largeobj] = max(szs);
largeobjidx = objs{largeobj};
largeobj = zeros(size(volthreshed)); largeobj(largeobjidx)=1;
voldist = bwdist(~largeobj);
figure; imshow(voldist(:,:,slc), []);
voldistcomp = -voldist; voldistcomp(~largeobj)=Inf;
Lwater = watershed(voldistcomp);
Lwater(~largeobj)=0;
figure; imshow(label2rgb(Lwater(:,:,slc)));
%%
Lwater2 = watershed(volhist);
max(Lwater2(:))