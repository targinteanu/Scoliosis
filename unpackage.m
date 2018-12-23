function vol = unpackage(packedvol, filter)

if nargin < 2
    filter = false;
end

sz = size(packedvol.Volume); 
rows = sz(1); cols = sz(2); ht = sz(3); % fix for 2D volume 

top = false(rows, cols, packedvol.minz-1);
bot = false(rows, cols, packedvol.maxz - ht - packedvol.minz + 1);

V = packedvol.Volume;

if filter
    V = medfilt3(V);
    V = imfill(V, 'holes');
end

vol = cat(3, top, V, bot);