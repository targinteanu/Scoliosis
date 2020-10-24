fp = uigetdir; 
[vol, imwidth, imheight, imdepth] = dicomreadvol(fp);

load([fp, '\vertebrae.mat']);

Vertebrae = Vertebrae(2:end); % if Vertebrae{1} is "none"

tf = false(size(Vertebrae));
for n = 1:length(Vertebrae)
    sz = size(Vertebrae{n}.Volume);
    tf(n) = length(sz) > 2;
end
vertebrae = Vertebrae(tf);

for n = 1:length(vertebrae)
    vertebrae{n} = unpackage(vertebrae{n}, true);
end

%%
figure; histogram(vol(:));

% patient 22
%%{
sliceToShow = 145;
slc = vol(:,:,sliceToShow);
figure; imshow(slc);
slc = slc - .24;
slc = slc.*(slc > 0);
seR = 5;
mx = .4;
%}
% patient 15
%{
sliceToShow = 32;
slc = vol(:,:,sliceToShow);
%slc = slc + .08;
slc = slc.*(slc > 0);
seR = 7;
%}
% patient 6
%{
sliceToShow = 30;
slc = vol(:,:,sliceToShow);
slc = slc - .24;
slc = slc.*(slc > 0);
seR = 6;
mx = .4-.24;
%}
% patient 28
%{
sliceToShow = 32;
slc = vol(:,:,sliceToShow);
slc = slc - .2;
slc = slc.*(slc > 0);
seR = 7;
mx = .6-.2;
%}

slc(slc > mx) = mx;

seg = zeros(imheight, imwidth);
for n = 1:length(vertebrae)
    v = vertebrae{n}(:,:,sliceToShow);
    seg(v(:)) = n;
end

figure; subplot(1,2,1); imshow(slc, []);
subplot(1,2,2); imshow(label2rgb(seg));

verts = unique(seg(:))

%% enlarge 
scl = 2;
slc_large = imresize(slc, scl, 'bilinear');
seg_large = imresize(seg,scl, 'nearest');

%% brighten
mapbrightness = @(I,k) 1-exp(-k*I);

figure; 
subplot(2,2,1); histogram(slc_large(:)); 
slc_bright = mapbrightness(slc_large,3); 
subplot(2,2,2); histogram(slc_bright(:));
subplot(2,2,3); imshow(slc_large, []);
subplot(2,2,4); imshow(slc_bright, []); 

%% sharpening matrix
figure; 
amt = 30:5:40;
rad = 1:.25:1.5;
for i = 1:length(amt)
    for j = 1:length(rad)
        spidx = (j-1)*length(rad) + i;
        subplot(length(rad), length(amt), spidx);
        
        slc_sharp = imsharpen(slc_bright, 'Amount', amt(i), 'Radius', rad(j));
        slc_sharp(seg_large==0) = slc_bright(seg_large==0);
        slc_sharp = slc_sharp.*(slc_sharp > 0);
        
        imshow(slc_sharp, []);
    end
end

%% sharpen 

%{
k = cell(1,6);
k{1} = [0,0,-1,0,0;0,-1,-2,-1,0;-1,-2,16,-2,-1;0,-1,-2,-1,0;0,0,-1,0,0];
k{2} = [0,0,-1,0,0;0,-1,-2,-1,0;-1,-2,17,-2,-1;0,-1,-2,-1,0;0,0,-1,0,0];
k{3} = [0,-1,0;-1,4,-1;0,-1,0];
k{4} = [0,-1,0;-1,5,-1;0,-1,0];
k{5} = [-1,-1,-1;-1,8,-1;-1,-1,-1];
k{6} = [-1,-1,-1;-1,9,-1;-1,-1,-1];
figure;
for i = 1:length(k)
    subplot(2,4,i);
    imshow(conv2(slc_large, k{i}), []);
end
%}

slc_sharp = imsharpen(slc_bright, 'Amount', 40, 'Radius', 1.25, 'Threshold', 0);
%slc_sharp = conv2(slc_large,k{4}); slc_sharp = slc_sharp(2:end-1,2:end-1);
slc_sharp(seg_large==0) = slc_bright(seg_large==0);
slc_sharp = slc_sharp.*(slc_sharp > 0);

figure; subplot(1,2,1); imshow(slc_bright, []); %subplot(2,2,3); histogram(slc(:));
subplot(1,2,2); imshow(slc_sharp, []); %subplot(2,2,4); histogram(slc_sharp(:));

%% show seg
colr = {'r', 'b', 'g', 'c', 'm', 'y'};
lsty = {'-', '--', ':', '-.'};
figure; imshow(slc_sharp, []); hold on;
se = strel('sphere', seR);
c = 1;
for n = 1:length(vertebrae)
    v = vertebrae{n}(:,:,sliceToShow);
    if sum(v(:))
        v = imdilate(v, se);
        v = imfill(v,'holes');
        v = imresize(v,scl, 'nearest');
        visboundaries(v, 'Color', colr{c}, 'LineStyle', lsty{c}, 'LineWidth', 1.5);
        c = c+1;
    end
end