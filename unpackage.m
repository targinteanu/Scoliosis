function vol = unpackage(packedvol, filter)
% reconstructs the original volume matrix that has been packaged by the
% package function. 
% packedvol: struct output by the function package. 
% filter: optional boolean input that specifies whether to filter the
%         volume matrix before unpackaging, which takes less time. If true,
%         the object will be smoothed with medfilt3 and holes will be
%         filled with imfill. If unspecified or false, no modifications
%         will be made to the volume before unpackaging. 

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