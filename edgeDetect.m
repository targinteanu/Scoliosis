fp = 'C:\Users\Toren\Desktop\scoliosis\Patient 17\Lumbar Spine C- Ct\L-SPINE - 2'; % pat 17 control
load([fp, '\background_import.mat']);

vol = dicomreadvol(fp); 

minval = min([vol(:); bg(:)]); maxval = max([vol(:); bg(:)]) - minval;
bg = bg - minval; bg = bg/maxval;
vol = vol - minval; vol = vol/maxval; % vol from 0 to 1

vol2 = vol - bg; 

volHE = histeq(vol); vol2HE = histeq(vol2);

%%
%figure; histogram(vol(:), 50); hold on; histogram(vol2, 50);
%volBW = (vol > .75); vol2BW = vol2 > .005; 
volBW = histThreshold(vol); vol2BW = histThreshold(vol2);
im1 = show3dBW(volBW, 100, [3,2,1]); im2 = show3dBW(volBW, 100, [3,1,2]);
im3 = show3dBW(vol2BW, 100, [3,2,1]); im4 = show3dBW(vol2BW, 100, [3,1,2]);
figure; imshow([im1, im2; im3, im4]);
%%
%{
figure; histogram(volHE(:), 50); hold on; histogram(vol2HE, 50);
volHEBW = (volHE > .96); vol2HEBW = vol2HE > .96; 
im1 = show3dBW(volHEBW, 100, [3,2,1]); im2 = show3dBW(volHEBW, 100, [3,1,2]);
im3 = show3dBW(vol2HEBW, 100, [3,2,1]); im4 = show3dBW(vol2HEBW, 100, [3,1,2]);
figure; imshow([im1, im2; im3, im4]);
%}

%%
%edgekernel(:,:,1) = -1*ones(3); edgekernel(:,:,3) = -1*ones(3); 
%    edgekernel(:,:,2) = [-1, -1, -1; -1, 26, -1; -1, -1, -1]; 
edgekernel(:,:,3) = [0,0,-1,0,0; 0,-1,-2,-1,0; -1,-2,30,-2,-1; 0,-1,-2,-1,0; 0,0,-1,0,0];
    edgekernel(:,:,2) = [0,0,0,0,0; 0,0,-1,0,0; 0,-1,-2,-1,0; 0,0,-1,0,0; 0,0,0,0,0];
    edgekernel(:,:,1) = [0,0,0,0,0,; 0,0,0,0,0; 0,0,-1,0,0; 0,0,0,0,0; 0,0,0,0,0];
    edgekernel(:,:,4) = edgekernel(:,:,2); edgekernel(:,:,5) = edgekernel(:,:,1);
edg = convn(vol, edgekernel); 
edg = edg-min(edg(:)); edg = edg/max(edg(:)); 
edg2 = histeq(edg);
%figure; histogram(edg(:), 100); hold on; histogram(edg2(:), 100);
%%
%edgBW = (edg > .44); 
%edgBW = (edg2 > .5);
edgBW = histThreshold(edg);
im3 = show3dBW(edgBW, 100, [3,2,1]); im4 = show3dBW(edgBW, 100, [3,1,2]);
figure; imshow([im3, im4]);

%%
radkernel = repmat(...
    [0,0,-1,0,0; 0,-1,-2,-1,0; -1,-2,80,-2,-1; 0,-1,-2,-1,0; 0,0,-1,0,0], ...
    [1,1,5]);
rad = convn(vol, radkernel); 
rad = rad-min(rad(:)); rad = rad/max(rad(:)); 
rad2 = histeq(rad);
%figure; histogram(rad(:), 100); hold on; histogram(rad2(:), 100);

%radBW = rad > .75;
radBW = histThreshold(rad);
im3 = show3dBW(radBW, 100, [3,2,1]); im4 = show3dBW(radBW, 100, [3,1,2]);
figure; imshow([im3, im4]);

%%
horkernel(:,:,1) = -1*ones(5); horkernel(:,:,2) = -2*ones(5); 
    horkernel(:,:,3) = 6*ones(5);
    horkernel(:,:,4) = horkernel(:,:,2); horkernel(:,:,5) = horkernel(:,:,1); 
hor = convn(vol, horkernel);
hor = hor-min(hor(:)); hor = hor/max(hor(:));
horBW = histThreshold(hor); 
im3 = show3dBW(horBW, 100, [3,2,1]); im4 = show3dBW(horBW, 100, [3,1,2]);
figure; imshow([im3, im4]);

figure; imshow([show3dBW(horBW, 50), show3dBW(radBW, 50)])

%%
%{
CC = bwconncomp(volBW); 
minsz = 10000; % pixels 
sz = arrayfun(@(i) length(CC.PixelIdxList{i}), 1:length(CC.PixelIdxList));
CC.PixelIdxList = CC.PixelIdxList(sz > minsz); 
CC.NumObjects = length(CC.PixelIdxList);
L = labelmatrix(CC); 

figure; 
for n = 1:CC.NumObjects
    subplot(CC.NumObjects,1,n);
    im3 = show3dBW(L==n, 100, [3,2,1]); im4 = show3dBW(L==n, 100, [3,1,2]);
    imshow([im3,im4]);
end
%}