n = 20; %c = 10; %optimal c between 10 and 11
cs = 140; 
i = 1;
imbalance = zeros(size(cs));
for c = cs

v = vertebrae{n}; 
[a1, a2, b1, b2] = PlaneOfSymmetry3(cm(n,:), sp(n,:));
c = c + 256;
a = a2-a1*c;
b = b2-b1*c;
pl = false(size(v));

for z = 1:199
    for y = 1:512
        x = a*z + b*y + c;
        x = floor(x);
        if (x > 0) & (x <= 512)
            pl(y, x, z) = true;
        end
    end
end

vpl = v | pl;

%%

yxz = find(v);
[y, x, z] = ind2sub(size(v), yxz);

idx = (x > (a*z + b*y + c));

xR = x(idx); yR = y(idx); zR = z(idx);
Rpts = (zR-1)*512*512 + (xR-1)*512 + yR;
%R = false(size(v)); R(Rpts) = true;

idx = (x < (a*z + b*y + c));
xL = x(idx); yL = y(idx); zL = z(idx);
Lpts = (zL-1)*512*512 + (xL-1)*512 + yL;
%L = false(size(v)); L(Lpts) = true;

imbalance(i) = length(Rpts) - length(Lpts);
i = i+1;
end