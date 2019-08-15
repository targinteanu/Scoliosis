function [imgthresh, hithresh, lothresh] = histThreshold(img, pNot, threshRatio)

if (nargin < 3)
    threshRatio = .4;
if (nargin < 2)
    pNot = .95;
end
end

numpix = 1; 
for s = size(img)
    numpix = numpix * s;
end
counts = imhist(img, 64);
hithresh = find(cumsum(counts) > pNot*numpix, 1, 'first') / 64;
lothresh = threshRatio * hithresh;
imgthresh = img > hithresh; 

end