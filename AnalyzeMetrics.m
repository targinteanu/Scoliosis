[num, txt] = xlsread('Writhe-pre-post_new-metrics_2-1.xlsx');
N = 32;
num = num(1:N, :);
XYZ = num(:, 13:63); 
ROT = num(:, 64:end);
IDCHOP = txt(2:(N+1),1); ID = txt(2:(N+1),3);

cluster_shape = num(:,1); cluster_writhe = num(:,2); 
    cluster_tor = num(:,3); cluster_twist = num(:,4); cluster_writhetwist = num(:,5);
    cluster = num(:,1:5);
    
writhe = num(:,6); abswrithe = num(:,7); 
tor1 = num(:,8); tor2 = num(:,9); torglob = num(:,10); 
twist = num(:,11); writhetwist = num(:,12);

%% correlations
%varnames2 = {'wr', 'tw'}; 
%varidx = [6,11];
%varnames2 = {'wr', 'tw', 'twr', 'tor', 'tor1'};
%varidx = [6, 11, 12, 10, 8];
varnames2 = {'wr', 'awr', '\tau_1', '\tau_2', '\tau_K', 'tw', 'twr'};
varidx = 6:12;

X = num(:,varidx); 
X(isnan(X))=0;
RHO = zeros(length(varnames2)); 
RHO(1,:) = arrayfun(@(x) corr(X(:,x), writhe), 1:length(varnames2)); 
RHO(2,:) = arrayfun(@(x) corr(X(:,x), twist), 1:length(varnames2)); 
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

%% correlation and p value of twist and writhe
[rho, p] = corr(writhe, twist); 
%[rho1m, p1m] = corr(writhe(cluster_shape==1), twist(cluster_shape==1));
%[rho2m, p2m] = corr(writhe(cluster_shape==2), twist(cluster_shape==2));
%[rho1w, p1w] = corr(writhe(cluster_writhe==1), twist(cluster_writhe==1));
C = [cluster_shape, cluster_writhe, cluster_twist];
RHO_P = cell(size(C,2), 2);
for r = 1:size(C,2)
    c = C(:,r);
    [rho1, p1] = corr(writhe(c==1), twist(c==1));
    [rho2, p2] = corr(writhe(c==2), twist(c==2));
    RHO_P{r,1} = [rho1; rho2]; RHO_P{r,2} = [p1; p2];
end
RHO_P = cell2mat(RHO_P);

%% mean, SD, p value
var = {'manual', 'writhe', 'torsion', 'twist', 'tw+wr'}; 

T = zeros(length(var), length(varnames2));
for v = 1:size(T,1)
    for j = 1:size(T,2)
        c = cluster(:,v);
        [~,T(v,j)] = ttest2(X(c==1,j), X(c==2,j));
    end
end

figure; heatmap(varnames2, var, T); 
title('p value'); xlabel('variable'); ylabel('clustering based on:');

AS = cell(length(var), length(varnames2));
for v = 1:size(AS,1)
    for j = 1:size(AS,2)
        c = cluster(:,v);
        avg = [mean(X(c==1,j)); mean(X(c==2,j))]; sd = [std(X(c==1,j)); std(X(c==2,j))];
        AS{v,j} = [avg, sd];
    end
end
AS = cell2mat(AS);

name1 = cell(1, 2*length(var)); 
for i = 1:length(var)
    i1 = (i-1)*2 + 1; i2 = (i-1)*2 + 2;
    name1{i1} = [var{i}, ' avg']; name1{i2} = [var{i}, ' std'];
end
name2 = cell(1, 2*length(varnames2));
for i = 1:length(varnames2)
    i1 = (i-1)*2 + 1; i2 = (i-1)*2 + 2;
    name2{i1} = [varnames2{i}, ' g1']; name2{i2} = [varnames2{i}, ' g2'];
end

figure; heatmap(name2, name1, AS); 
title('avg/std'); xlabel('variable'); ylabel('clustering based on:');

%% clustering 
var1 = twist; 
cluster_var1 = kmeans(var1, 2); 
acc = cluster_shape == cluster_var1; acc = sum(acc)/length(acc); 
acc
var2 = writhe; 
cluster_var2 = kmeans(var2, 2); 
acc = cluster_shape == cluster_var2; acc = sum(acc)/length(acc); 
acc
cluster_var12 = kmeans([var1 var2], 2); 
acc = cluster_shape == cluster_var12; acc = sum(acc)/length(acc); 
acc

%% groups 1 & 2 averages 
c = cluster_shape;
numclusters = 2;
CM = cell(1,numclusters); rot = cell(1,numclusters);
%%{
for i = 1:numclusters
    % for each cluster
    xyz = mean(XYZ(c==i,:));
    cm = arrayfun(@(j) xyz(j:3:end)', 1:3, 'UniformOutput', false);
    CM{i} = cell2mat(cm);
    rot{i} = mean(ROT(c==i,:))';
end
IDS = arrayfun(@(i) ['Avg group ',num2str(i)], 1:numclusters, 'UniformOutput', false);
%}
%{
ex = [17, 30];
for i = 1:numclusters
    xyz = XYZ(ex(i),:); 
    cm = arrayfun(@(j) xyz(j:3:end)', 1:3, 'UniformOutput', false);
    CM{i} = cell2mat(cm);
    rot{i} = ROT(ex(1),:)';
end
IDS = arrayfun(@(i) [ID{i},' (group ',num2str(c(i)),')'], ex, 'UniformOutput', false);
%}

vlabels_startend = cell(size(rot{1})); vlabels_startend{1} = 'T1'; vlabels_startend{end} = 'L5';
vlabels = [arrayfun(@(v) ['T' num2str(v)], 1:12, 'UniformOutput', false), ...
    arrayfun(@(v) ['L' num2str(v)], 1:5, 'UniformOutput', false)];
lbl = {'b) ', 'c) ', 'd) ', 'e) ', 'f) ', 'g) '};

figure('Position', [0 50 1500 500]); 
subplot(1,numclusters+1,1); 
LS = {{':b', 'ob', ':ob'}, {'--k','sk','--sk'}}; 
DS = {{':m', 'om'},        {'--r','sr'}};
for i = 1:numclusters
    %plot3dSpine(CM{i}, rot{i}, LS{i}, DS{i}); hold on; grid on;
    plot3(CM{i}(:,1), CM{i}(:,2), CM{i}(:,3), LS{i}{3}, 'LineWidth', 2); 
    hold on; grid on;
    text(CM{i}(:,1) - 20, CM{i}(:,2), CM{i}(:,3), vlabels_startend);
    title('a) Axial View', 'FontSize', 13);
end
view(2); 
xlabel('x'); ylabel('y');
legend(IDS, 'Location', 'southoutside');

for i = 1:numclusters
    subplot(1,numclusters+1,i+1); 
    plot3dSpine(CM{i}, rot{i}, LS{i}, DS{i});
    hold on; text(CM{i}(:,1), CM{i}(:,2) + 30, CM{i}(:,3), vlabels);
    grid on; view([-70,8]);
    xlabel('x'); ylabel('y'); zlabel('z');
    Tw = getTwist(CM{i}, rot{i}); Tw = round(Tw, 3);
    Wr = levittWrithe(CM{i}); Wr = round(Wr, 3);
    title([lbl{i}, 'Twist = ' num2str(Tw) ', Writhe = ' num2str(Wr)], 'FontSize', 13);
end