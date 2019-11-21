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

%% align coordinate systems 
clear ax;
% get ref pts 
figure; ax(1)=subplot(1,2,1); imshow(cor); ax(2)=subplot(1,2,2); imshow(sag);
[x,z1] = getpts(ax(1)); [y,z2] = getpts(ax(2));
% ignore excess 
zL = min(length(z1), length(z2));
x = x(1:zL); y = y(1:zL); z1 = z1(1:zL); z2 = z2(1:zL);
% find best fit 
c = [z1, ones(size(z1))]\z2; % z2 = c(1)*z1 + c(2)
%c = int16(c);

%% polynomial & spline fitting - works well
figure; 
polyfun = @(p0, p1, p2, p3, p4, p5, x) ...
    p0 + p1*x + p2*x.^2 + p3*x.^3 + p4*x.^4 + p5*x.^5; % + ...
%    p6*x.^6 + p7*x.^7 + p8*x.^8 + p9*x.^9 + p10*x.^10;
npoly = 10;
mm = 40;
polyval = @(p, x) sum( p.*(x'.^fliplr((1:length(p))-1)), 2);

[rsag, csag] = find(sagspine);
[fosag, gofsag] = fit(rsag, csag, 'smoothingspline', 'Weights', sagspine(find(sagspine(:)))); 
%[fpsag, gofsag] = fit(rsag, csag, polyfun, 'Weights', sagspine(find(sagspine(:))));
polysag = polyfit(rsag, csag, npoly);
img = zeros([size(sagspine),3]); img(:,:,1) = sag2; 
img(:,:,3) = sag2-3*sagspine; img(:,:,2) = img(:,:,3);
subplot(1,3,1); imshow(img); hold on; title('sag');
%plot(feval(fpsag,min(rsag):max(rsag)), min(rsag):max(rsag));
xsag = min(rsag):max(rsag); clsag = ppval(xsag, fosag.p);
clsag = movmean(clsag, mm); 
plot(clsag, xsag, 'b');
%plot(polyval(polysag, xsag), xsag);

[rcor, ccor] = find(corspine);
[focor, gofcor] = fit(rcor, ccor, 'smoothingspline', 'Weights', corspine(find(corspine(:))));
%[fpcor, gofcor] = fit(rcor, ccor, polyfun, 'Weights', corspine(find(corspine(:))));
polycor = polyfit(rcor, ccor, npoly);
img = zeros([size(corspine),3]); img(:,:,1) = cor2; 
img(:,:,3) = cor2-3*corspine; img(:,:,2) = img(:,:,3);
subplot(1,3,2); imshow(img); hold on; title('cor');
%plot(feval(fpcor,min(rcor):max(rcor)), min(rcor):max(rcor));
xcor = min(rcor):max(rcor); clcor = ppval(xcor, focor.p);
clcor = movmean(clcor, mm);
plot(clcor, xcor, 'b');
%plot(polyval(polycor, xcor), xcor);

% spline fit with aligned z
% z2 = c(1)*z1 + c(2)
% z1 = cor; z2 = sag;

[rsag, csag] = find(sagspine);
[rcor, ccor] = find(corspine);
rcor = c(1)*rcor + c(2); % rcor now aligned to rsag
z = max(min(rsag), min(rcor)):min(max(rsag), max(rcor));

[fosag, gofsag] = fit(rsag, csag, 'smoothingspline', 'Weights', sagspine(find(sagspine(:)))); 
clsag = ppval(z, fosag.p);
clsag = movmean(clsag, mm); 

[focor, gofcor] = fit(rcor, ccor, 'smoothingspline', 'Weights', corspine(find(corspine(:))));
clcor = ppval(z, focor.p);
clcor = movmean(clcor, mm);

subplot(1,3,3); plot3(clcor, clsag, z, 'b'); grid on; 
title(['Writhe = ' num2str(levittWrithe([clcor; clsag; z]'))]);

%% display results of coordinate system alignment - does not work
%{
corscl = cor2; sagscl = sag2;
z1scl = z1; z2scl = z2;
if c(2) < 0
    % cor higher than sag 
    sagscl = [zeros(-int16(c(2)), size(sag,2)); sag];
    z2scl = z2scl - c(2);
elseif c(2) > 0
    % sag higher than cor 
    corscl = [zeros(int16(c(2)), size(cor,2)); cor];
    z1scl = z1scl + c(2);
end
if c(1) > 1
    % stretch cor 
    corscl0 = corscl(1:int16((size(sag,1)/c(1))), :);
    corscl = zeros(size(sagscl));
    for col = 1:size(corscl,2)
        corscl(:,col) = interp(corscl0(:,col), c(1)); 
    end
    z1scl = z1scl * c(1);
elseif c(1) < 1
    % stretch sag
    sagscl0 = sagscl(1:int16((size(cor,1)*c(1))), :);
    sagscl = zeros(size(corscl));
    for col = 1:size(sagscl,2)
        sagscl(:,col) = interp(sagscl0(:,col), 1/c(1)); 
    end
    z2scl = z2scl / c(1);
end
figure; 
subplot(1,2,1); imshow(corscl); hold on; plot(x, z1scl, '*');
subplot(1,2,2); imshow(sagscl); hold on; plot(y, z2scl, '*');
%}

%% segment 
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

%% conv with different kernels 
img = sag2; 

kernels = {...
    [-1 -1 -1; -1 9 -1; -1 -1 -1],...
    [-1 2 1], [1 2 -1], [1; 2; -1], [-1; 2; 1], ...
    [-1 3 -1], [-1; 3; -1], ...
    [0 .25 0; .25 1 .25; 0 .25 0]};
imf = cell(size(kernels));
for i = 1:length(kernels)
    k = kernels{i};
    imf{i} = imfilter(img, k, 'conv', 'replicate');
end
figure; imshow(cell2mat(imf));

%% use kmeans on pixels 

R = repmat((1:size(img,1))', [1,size(img,2)]); C = repmat((1:size(img,2)), [size(img,1),1]);
R = R/max(R(:)); C = C/max(C(:));

G = 1:size(img,2); G = (1./sqrt(2*pi*std(G).^2)) .* exp(-((G-mean(G)).^2)./(2*std(G).^2));
G = repmat(G, [size(img,1),1]);

%X = [img(:), R(:), C(:)];
X = [img(:), G(:), imf{1}(:), imf{6}(:), imf{7}(:), imf{8}(:)];
idx = kmeans(X, 7); 
imgout = zeros(size(img)); imgout(:) = idx;
figure; imshow(label2rgb(imgout));
%figure; imshow([img, imgout/max(imgout(:))])

%% doesn't work well 
%{
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

%
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
%}