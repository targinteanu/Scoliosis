function packed = package(vol)
% Packages a large volume matrix into a struct that uses less memory. 
% vol is a volume matrix that has empty regions on the top and/or bottom. 
% packed.Volume is the region of vol with nonzero information. 
% packed.minz indicates at what height in vol packed.Volume begins. 
% packed.maxz indicates where vol ends.  

packed = struct;
maxz = size(vol); packed.maxz = maxz(3);
zidx = find(squeeze(sum(sum(vol))));
if isempty(zidx)
    packed.minz = 1;
    endz = 1;
else
    packed.minz = zidx(1); 
    endz = zidx(end);
end
packed.Volume = vol(:,:, packed.minz:endz);

end