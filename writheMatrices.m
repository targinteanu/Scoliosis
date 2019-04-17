idxs = 1:33; 
numcurves = zeros(size(idxs));
for idx = idxs
    writheM
    cc = bwconncomp(M > 0); 
    numP = cc.NumObjects;
    cc = bwconncomp(M < 0); 
    numN = cc.NumObjects;
    numcurves(idx) = min(numN, numP)/2;
end
figure; plot(numcurves)