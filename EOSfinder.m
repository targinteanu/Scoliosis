%[fn, fp] = uigetfile;
%img = dicomread([fp '\' fn]);
%ifo = dicominfo([fp '\' fn])
%figure; imshow(img, [])
fp = uigetdir; 
fs = 1:4;
for f = fs
   img = dicomread([fp '\(' num2str(f) ')']); 
   sz = size(img);
       figure; imshow(img, []);
       title([num2str(f) ' ' num2str(sz(1)) ' ' num2str(sz(2))]);
end