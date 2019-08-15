function imtoshow = show3dBW(bwvol, numLevels, perm)
% displays a 3D projection of BW object 
%   numLevels: number of slices of the object to show 
%   perm: permute the object before showing 

if (nargin > 2)
    bwvol = permute(bwvol, perm);
elseif (nargin < 2)
    numLevels = 10;
end

sz = size(bwvol); imdepth = sz(3);

cm = CenterOfMass3(bwvol); 
k1 = min(imdepth, floor(cm(3) + numLevels/2));
imtoshow = (1/numLevels)*bwvol(:,:,k1);
for k = 1:(numLevels-1)
    imtoshow(bwvol(:,:,max(1,k1-k))) = 1/numLevels + k/numLevels;
end

end