%[fn, fp] = uigetfile;
%img = dicomread([fp '\' fn]);
%ifo = dicominfo([fp '\' fn])
%figure; imshow(img, [])
fp = uigetdir; 
f = 1;
while exist([fp '\(' num2str(f) ')'], 'file')
   img = dicomread([fp '\(' num2str(f) ')']); 
   sz = size(img);
       figure; imshow(img, []);
       title([num2str(f) ' ' num2str(sz(1)) ' ' num2str(sz(2))]);
   f = f+1;
end