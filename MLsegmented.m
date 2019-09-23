load('vertebrae_pat6.mat'); 

tf = false(size(Vertebrae));
for n = 1:length(Vertebrae)
    sz = size(Vertebrae{n}.Volume);
    tf(n) = length(sz) > 2;
end
vertebrae = Vertebrae(tf);

for n = 1:length(vertebrae)
    vertebrae{n} = unpackage(vertebrae{n}, true);
end

%%

v1 = vertebrae{4};
v2 = vertebrae{5};
v = 2*v1 + v2;

k1 = find(v1); 
s1 = size(v1); 
z1 = floor(k1 / (s1(1)*s1(2)));
k1 = k1 - (s1(1)*s1(2)*z1);
y1 = floor(k1 / s1(1));
x1 = k1 - s1(1)*y1;
x1 = x1+1; y1 = y1+1; z1 = z1+1;

k2 = find(v2); 
s2 = size(v2); 
z2 = floor(k2 / (s2(1)*s2(2)));
k2 = k2 - (s2(1)*s2(2)*z2);
y2 = floor(k2 / s2(1));
x2 = k2 - s2(1)*y2;
x2 = x2+1; y2 = y2+1; z2 = z2+1;

spcTrn = [x1, y1, z1; x2, y2, z2]; vertTrn = [1*ones(size(x1)); 2*ones(size(x2))];

%%
olp = sum(sum(v1)) & sum(sum(v2)); olp = squeeze(olp); olp = find(olp);

[spaceX, spaceY, spaceZ] = meshgrid(1:(s1(1)), 1:(s1(2)), min(olp):max(olp));
spcTest = [spaceX(:), spaceY(:), spaceZ(:)];

vertML = knnclassify(spcTest, spcTrn, vertTrn, 5);

%%
vML = reshape(vertML, s1(1), s1(2), length(olp)); vML = permute(vML, [2 1 3]);
X = zeros(size(v)); X(:,:,olp) = vML; 
vv = v+X;