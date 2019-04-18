idxs = 1:33; 
numcurves = zeros(size(idxs));
for idx = idxs
        xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    
cm = [x; y; z]';

    writheM
    cc = bwconncomp(M > 0); 
    numP = cc.NumObjects;
    cc = bwconncomp(M < 0); 
    numN = cc.NumObjects;
    numcurves(idx) = min(numN, numP)/2;
end
figure; plot(numcurves)