function packed = package(vol)

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