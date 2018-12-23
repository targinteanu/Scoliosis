function [V, imwidth, imheight, imdepth] = dicomreadvol(filepath)

filebase = '\IM-0001-'; fileend = '-0001.dcm';
i = 1; 
filename = [filepath, filebase, num2str(i, '%04u'), fileend];
slice = dicomread(filename);
imwidth = size(slice); imheight = imwidth(1); imwidth = imwidth(2);
V = zeros(imheight, imwidth, 1000);
%vol = int16(vol);
maxvol = 4095;

%

while exist(filename, 'file')
    slice = dicomread(filename);
    V(:,:,i) = slice; 
    
    i = i+1;
    filename = [filepath, filebase, num2str(i, '%04u'), fileend];
end

imdepth = i;
V = V(:,:,1:i);
V = V/maxvol;

end