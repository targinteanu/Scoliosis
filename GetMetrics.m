[num, txt] = xlsread('Writhe-pre-post_new-metrics_1-25.xlsx');
N = 32;
num = num(1:N, :);
XYZ = num(:, 13:63); 
ROT = num(:, 64:end);

cluster_shape = num(:,1); cluster_writhe = num(:,2); 
    cluster_tor = num(:,3); cluster_twist = num(:,4); cluster_writhetwist = num(:,5);
    
writhe = num(:,6); abswrithe = num(:,7); 
tor1 = num(:,8); tor2 = num(:,9); torglob = num(:,10); 
twist = num(:,11); writhetwist = num(:,12);

%% example plot
pat = 3;

x = XYZ(pat, 1:3:end)'; 
y = XYZ(pat, 2:3:end)';
z = XYZ(pat, 3:3:end)';
theta = ROT(pat, :)'; %theta = -theta;

figure; plot3dSpine([x,y,z], theta);

%% re-evaluate metrics 
newwrithe = arrayfun(@(p) ...
    levittWrithe([XYZ(p,1:3:51); XYZ(p,2:3:51); XYZ(p,3:3:51)]'), 1:N);
newtwist = arrayfun(@(p) ...
    getTwist([XYZ(p,1:3:51); XYZ(p,2:3:51); XYZ(p,3:3:51)]', ROT(p,:)'), 1:N);
newtorsion = arrayfun(@(p) ...
    kadouryTorsion([XYZ(p,1:3:51); XYZ(p,2:3:51); XYZ(p,3:3:51)]'), 1:N);

newwrithe = newwrithe'; newtwist = newtwist'; newtorsion = newtorsion';
newwrithetwist = newtwist + newwrithe;
newtorsion(isnan(newtorsion)) = 0; torglob(isnan(torglob)) = 0;

torsionchange = sum(abs(newtorsion - torglob))
twistchange = sum(abs(newtwist - twist))
writhechange = sum(abs(newwrithe - writhe))
writhetwistchange = sum(abs(writhetwist - newwrithetwist))

%% clustering 
cluster_newtor = kmeans(newtorsion, 2);
cluster_newtwist = kmeans(newtwist, 2); 
cluster_newwrithe = kmeans(newwrithe, 2);
cluster_newwrithetwist = kmeans(newwrithetwist, 2);

acc = @(clus) mean(clus == cluster_shape);