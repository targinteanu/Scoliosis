cor = imread('eso001cor.png'); sag = imread('eso001sag.png'); 
cor = rgb2gray(cor); sag = rgb2gray(sag);
%figure; corroi = roipoly(cor); 
%figure; sagroi = roipoly(sag); 
load('esotestroi.mat');

%%
cortrim = cor.*uint8(corroi); sagtrim = sag.*uint8(sagroi); 
figure; imshow(cortrim); figure; imshow(sagtrim);

%%
sag2 = double(sag)/255; cor2 = double(cor)/255;
%{
im = showreconstruct(sag2); sagbone = sag2-im; 
im = showreconstruct(cor2); corbone = cor2-im; 
%}
sagim = imreconstruct(double(~sagroi), sag2); sagbone = sag2-sagim; 
corim = imreconstruct(double(~corroi), cor2); corbone = cor2-corim; 
sagbone = adapthisteq(sagbone); corbone = adapthisteq(corbone); 

%sagspine = sagbone; corspine = corbone;
sagspine = sagbone.*sagroi; corspine = corbone.*corroi; 
spine = [sagspine corspine]; 
spinethresh = [spine > .03*max(spine(:)), spine > .05*max(spine(:)), spine > .07*max(spine(:))];
%figure; imshow(spinethresh);
figure; imshow([sag2 sagspine cor2 corspine])
%}

%%
seg = corspine;
CC = bwconncomp(seg > .04);
minsz = 5; % pixels 
sz = arrayfun(@(i) length(CC.PixelIdxList{i}), 1:length(CC.PixelIdxList));
CC.PixelIdxList = CC.PixelIdxList(sz > minsz); 
CC.NumObjects = length(CC.PixelIdxList);
L = labelmatrix(CC); %figure; imshow(label2rgb(L));
CMs = zeros([size(seg), 3]); 
for pixIdx = CC.PixelIdxList
    obj = zeros(size(seg)); obj(pixIdx{:}) = 1;
    cm = CenterOfMass(obj); cm = uint16(cm); 
    CMs(cm(1), cm(2), :) = 1;
end

seg = sagspine;
CC = bwconncomp(seg > .04);
minsz = 5; % pixels 
sz = arrayfun(@(i) length(CC.PixelIdxList{i}), 1:length(CC.PixelIdxList));
CC.PixelIdxList = CC.PixelIdxList(sz > minsz); 
CC.NumObjects = length(CC.PixelIdxList);
Lsag = labelmatrix(CC); 
CMs_sag = zeros([size(seg), 3]); 
for pixIdx = CC.PixelIdxList
    obj = zeros(size(seg)); obj(pixIdx{:}) = 1;
    cm = CenterOfMass(obj); cm = uint16(cm); 
    CMs_sag(cm(1), cm(2), :) = 1;
end

figure; imshow([label2rgb(L), 255*CMs, label2rgb(Lsag), 255*CMs_sag]);

%%
centlinesSag = {double(~~Lsag), sagspine, CMs_sag(:,:,1)};
for i = 1:length(centlinesSag)
    img = centlinesSag{i};
    imgout = zeros(size(img));
    for r = 1:size(img,1)
        cm = CenterOfMass1(img(r,:));
        cm = uint16(cm);
        if cm
            imgout(r, cm) = 1;
        end
    end
    centlinesSag{i} = imgout;
end
figure; imshow(cell2mat(centlinesSag));

%%
orig = sagspine;
roughline = CenterOfMass1(orig);
rang = max(roughline)-min(roughline);

[pk1, loc1, pw1, pp1] = findpeaks(roughline, 'MinPeakProminence', rang/10);
[pk2, loc2, pw2, pp2] = findpeaks(-roughline, 'MinPeakProminence', rang/10); pk2=-pk2;

pk = [pk1; pk2]; loc = [loc1; loc2]; pw = [pw1; pw2]; pp = [pp1; pp2];
[loc,idx] = sort(loc); pk=pk(idx); pw=pw(idx); pp=pp(idx);

xq = 1:size(orig,1);
rlspl = spline(loc, pk, xq);

figure; imshow(orig); hold on; 
plot(roughline, xq, '.y');
plot(rlspl, xq, 'y');
errorbar(pk,loc, pw,pw, pp,pp, 'oy');