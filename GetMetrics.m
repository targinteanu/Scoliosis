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

%%
pat = 10;

x = XYZ(pat, 1:3:end)'; 
y = XYZ(pat, 2:3:end)';
z = XYZ(pat, 3:3:end)';
theta = ROT(pat, :)'; %theta = -theta;

figure; plot3dSpine([x,y,z], theta);

%%
checkwrithe = arrayfun(@(p) ...
    levittWrithe([XYZ(p,1:3:51); XYZ(p,2:3:51); XYZ(p,3:3:51)]'), 1:N);
checktwist = arrayfun(@(p) ...
    getTwist([XYZ(p,1:3:51); XYZ(p,2:3:51); XYZ(p,3:3:51)]', ROT(p,:)'), 1:N);
checkwrithe = checkwrithe'; checktwist = checktwist';
sum(abs(checktwist - twist))
sum(abs(checkwrithe - writhe))

%%
varnames2 = {'wr', 'tw'}; 
varidx = [6,11];
X = num(:,varidx); 
X(isnan(X))=0;
RHO = zeros(length(varnames2)); 
RHO(1,:) = arrayfun(@(x) corr(X(:,x), writhe), 1:length(varnames2)); 
RHO(2,:) = arrayfun(@(x) corr(X(:,x), checktwist), 1:length(varnames2)); 
cmp = [linspace(1, 1), linspace(1, 0); 
        linspace(0, 1), linspace(1, 1);
        linspace(1, 1), linspace(1, 1)]';
figure; heatmap(varnames2, varnames2, RHO, 'ColorLimits', [-1 1], 'Colormap', cmp); 
title('correlation between variables');

for group = 1:2
Xg = num(cluster_shape==group,varidx); 
Xg(isnan(Xg))=0;
RHO = zeros(length(varnames2)); 
for r = 1:length(varnames2)
    RHO(r,:) = arrayfun(@(x) corr(Xg(:,x), Xg(:,r)), 1:length(varnames2)); 
end

cmp = [linspace(1, 1), linspace(1, 0); 
        linspace(0, 1), linspace(1, 1);
        linspace(1, 1), linspace(1, 1)]';
figure; heatmap(varnames2, varnames2, RHO, 'ColorLimits', [-1 1], 'Colormap', cmp); 
title(['correlation between variables group ' num2str(group)]);
end