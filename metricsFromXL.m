[num, txt] = xlsread('Writhe-pre-post_new-metrics_10-14.xlsx');
N = 32;
num = num(1:N, :);
XYZ = num(:, 13:end); 

cluster_shape = num(:,1); cluster_writhe = num(:,2); 
    cluster_tor = num(:,3); cluster_twist = num(:,4); cluster_writhetwist = num(:,5);
%writhe_preop = num(:,4); 
%writhe_postop = num(:,5);
%nonsurg = isnan(writhe_postop);

writhe = num(:,6); abswrithe = num(:,7); 
tor1 = num(:,8); tor2 = num(:,9); torglob = num(:,10); 
twist = num(:,11); writhetwist = num(:,12);

%%
var = {'manual', 'writhe', 'torsion', 'twist'}; 
cluster = num(:,1:4);
sp = 220; figure; 
for v = 1:length(var)
    %ax(v) = subplot(sp + v);
    subplot(sp + v);
    c = cluster(:,v);
    %x = cell(1,2); y = cell(1,2); z = cell(1,2);
    lgd = cell(1,2);
    for group = 1:2
        x = XYZ(c == group, 1:3:51);
        y = XYZ(c == group, 2:3:51);
        z = XYZ(c == group, 3:3:51);
        
        p = [mean(x); mean(y); mean(z)]';
        [tau, neutral, apical] = kadouryTorsion(p); taupts = [neutral, apical];
        Wr = levittWrithe(p); 
        
        lgd{group} = ['group ' num2str(group) ' | Wr=' num2str(Wr) ' | \tau=' num2str(tau)];
        
        %errorbar(mean(x), mean(y), std(y),std(y), std(x),std(x), '-o'); 
        plot3(mean(x), mean(y), mean(z), '-o'); 
        hold on;
        plot3(p(taupts,1), p(taupts,2), p(taupts,3), '*');
    end
    grid on; view([0,90]);
    lgd = {lgd{1}, 'neutral/apical', lgd{2}, 'neutral/apical'}; 
    legend(lgd);
    title(var{v}); 
end
%linkaxes(ax);

%% check values 
rotv = @(theta) [cos(theta), sin(theta), zeros(size(theta))];

checkwrithe = arrayfun(@(p) ...
    levittWrithe([XYZ(p,1:3:51); XYZ(p,2:3:51); XYZ(p,3:3:51)]'), 1:N);
checktwist = arrayfun(@(p) ...
    getTwist([XYZ(p,1:3:51); XYZ(p,2:3:51); XYZ(p,3:3:51)]', rotv(XYZ(p,52:end)')), 1:N);

sum(abs(checktwist' - twist))
sum(abs(checkwrithe' - writhe))

%%
varnames2 = {'wr', 'awr', '\tau_1', '\tau_2', '\tau_K', 'tw', 'twr'};
X = num(:,6:12); 
X(isnan(X))=0;
RHO = zeros(length(varnames2)); 
RHO(1,:) = arrayfun(@(x) corr(X(:,x), writhe), 1:7); 
RHO(2,:) = arrayfun(@(x) corr(X(:,x), abswrithe), 1:7); 
RHO(3,:) = arrayfun(@(x) corr(X(:,x), tor1), 1:7); 
RHO(4,:) = arrayfun(@(x) corr(X(:,x), tor2), 1:7); 
RHO(5,:) = arrayfun(@(x) corr(X(:,x), X(:,5)), 1:7); % torglob
RHO(6,:) = arrayfun(@(x) corr(X(:,x), twist), 1:7); 
RHO(7,:) = arrayfun(@(x) corr(X(:,x), writhetwist), 1:7); 
%RHO
%cmp = [1 0 0; 1 1 1; 0 0 1];
cmp = [linspace(1, 1), linspace(1, 0); 
        linspace(0, 1), linspace(1, 1);
        linspace(1, 1), linspace(1, 1)]';
figure; heatmap(varnames2, varnames2, RHO, 'ColorLimits', [-1 1], 'Colormap', cmp); 
title('correlation between variables');

%%
var = {'manual', 'writhe', 'torsion', 'twist', 'tw+wr'}; 
cluster = num(:,1:5);

T = zeros(length(var), length(varnames2));
for v = 1:size(T,1)
    for j = 1:size(T,2)
        c = cluster(:,v);
        [~,T(v,j)] = ttest2(X(c==1,j), X(c==2,j));
    end
end

figure; heatmap(varnames2, var, T); 
title('p value'); xlabel('variable'); ylabel('clustering based on:');

%%
varnames2 = {'wr', 'awr', '\tau_1', '\tau_2', '\tau_K', 'tw', 'twr'};
group = 1;
Xg = num(cluster_shape==group,6:12); 
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

%%
%{
pat = zeros(4, length(varnames2));
for r = 1:4
    pat(r, :) = arrayfun(@(v) mean(X(cluster(:,r)==1,v)) > mean(X(cluster(:,r)==2,v)), ...
        1:length(varnames2));
end
figure; heatmap(varnames2([1,6,5,3,4]), var, pat(:, [1,6,5,3,4]));
%}

%%

%{
writhes = zeros(N,1); abswrithes = zeros(N,1); 
decompwrithes = zeros(N,1); crosswrithe = zeros(N,1); crosswrithes = zeros(N,2);
topbottomwrithes = zeros(N,1);
torsions = zeros(N,1); torsionlocs = zeros(N,1);
torsions2 = zeros(N,1); torsionlocs2 = zeros(N,1);
twists = zeros(N,1); 
for idx = 1:N
    % get center points 
    x = XYZ(idx, 1:3:51); 
    y = XYZ(idx, 2:3:51); 
    z = XYZ(idx, 3:3:51); 
    cm = [x;y;z]';
    
    % get rotation vectors 
    theta = XYZ(idx, 52:end)'; 
    rotvector = [cos(theta), sin(theta), zeros(size(theta))];
        % this will be made orthogonal to dr in the twist function. 
    
    % get writhe 
    writhes(idx) = levittWrithe(cm);
    abswrithes(idx) = levittWritheAbs(cm);
%    decompwrithes(idx) = decompWrithe(cm); 
    [crosswrithe(idx), crosswrithes(idx,:)] = crossWrithe(cm);
    topbottomwrithes(idx) = levittWrithe(cm, [1 2 16 17]);
    
    % get twist
    twists(idx) = getTwist(cm, rotvector);
    
    % get torsion 
    q = 4; % 2 vertebrae above, 2 vertebrae below, current vertebra -> 5 points to fit each cubic
    vertebrae = (1+q):(length(x)-q);
    tau = arrayfun(@(v) lewinerTorsion(cm, v, q), vertebrae);
    [~,tauvert] = max(abs(tau)); torsions(idx) = tau(tauvert); % maximum torsion
    torsionlocs(idx) = tauvert+q;
    
    % get second-largest torsion
    tau2 = tau([1:(tauvert-1), (tauvert+1):end]);
    [~,tauvert2] = max(abs(tau2)); torsions2(idx) = tau2(tauvert2); 
    torsionlocs2(idx) = tauvert2 + q + (tauvert2>=tauvert);
    
%    idx/N
end
%}

%abstorsions = abs(torsions); abstorsions2 = abs(torsions2);

%{
%% hypothesis test
[h,p] = ttest(abs(writhes(~nonsurg)), abs(writhe_postop(~nonsurg)), 'Tail', 'right')
figure; boxplot(abs([writhes(~nonsurg), writhe_postop(~nonsurg)]), ...
    'Labels', {'pre-operative', 'post-operative'});
ylabel('Writhe Magnitude');

%%
% make sure writhes from table are the same as writhes recalculated 
% figure; plot(writhe_preop); hold on; grid on; plot(writhes);

% compare writhes pre and post op
figure; plot(writhes(~nonsurg), '.r'); hold on; grid on; plot(writhe_postop(~nonsurg), '.b');

% compare spines that received surgery to those that didn't 
figure; plot(find(nonsurg), abs(writhes(nonsurg)), 'xb'); hold on; grid on; 
plot(find(nonsurg), abswrithes(nonsurg), 'ob'); 
plot(find(nonsurg), abs(torsions(nonsurg)), '^b'); 
plot(find(~nonsurg), abs(writhes(~nonsurg)), 'xr'); 
plot(find(~nonsurg), abswrithes(~nonsurg), 'or'); 
plot(find(~nonsurg), abs(torsions(~nonsurg)), '^r'); 

% compare clusters 1 and 2 
figure; plot(find(shapecluster == 1), abs(writhes(shapecluster == 1)), 'xb'); hold on; grid on; 
plot(find(shapecluster == 1), abswrithes(shapecluster == 1), 'ob'); 
plot(find(shapecluster == 1), abs(torsions(shapecluster == 1)), '^b'); 
plot(find(shapecluster == 2), abs(writhes(shapecluster == 2)), 'xr'); 
plot(find(shapecluster == 2), abswrithes(shapecluster == 2), 'or'); 
plot(find(shapecluster == 2), abs(torsions(shapecluster == 2)), '^r');

%% 3d cluster plots
figure; 
plot3(writhes(shapecluster == 1), abswrithes(shapecluster == 1), ...
    abstorsions(shapecluster == 1), 'ob'); 
grid on; hold on; 
plot3(writhes(shapecluster == 2), abswrithes(shapecluster == 2), ...
    abstorsions(shapecluster == 2), '^r'); 
xlabel('writhe'); ylabel('absolute writhe'); zlabel('torsion');

%% 2d cluster plots
figure; 
subplot(2,2,1); 
plot(writhes(shapecluster == 1), abstorsions(shapecluster == 1), 'ob'); 
grid on; hold on; 
plot(writhes(shapecluster == 2), abstorsions(shapecluster == 2), '^r'); 
legend('Group 1', 'Group 2'); title('A) Writhe and Torsion'); 
xlabel('Writhe'); ylabel('Torsion Magnitude'); 

%{
subplot(2,2,2); 
plot(abs(writhes(shapecluster == 1)), abswrithes(shapecluster == 1), 'ob'); 
grid on; hold on; 
plot(abs(writhes(shapecluster == 2)), abswrithes(shapecluster == 2), '^r'); 
legend('Group 1', 'Group 2'); title('B) Writhe Magnitude and Absolute Writhe'); 
xlabel('Writhe Magnitude'); ylabel('Absolute Writhe'); 
%}
subplot(2,2,2); 
plot(abswrithes(shapecluster == 1), abstorsions(shapecluster == 1), 'ob'); 
grid on; hold on; 
plot(abswrithes(shapecluster == 2), abstorsions(shapecluster == 2), '^r'); 
legend('Group 1', 'Group 2'); title('B) Absolute Writhe and Torsion'); 
xlabel('Absolute Writhe'); ylabel('Torsion Magnitude'); 

% PCA clustering 
[C,S,~,~,P] = pca(XYZ(:,1:51));
S1 = S(find(shapecluster==1),:); S2 = S(find(shapecluster==2),:);
subplot(2,2,3); 
plot(S1(:,1), S1(:,2), 'ob'); grid on; hold on; 
plot(S2(:,1), S2(:,2), '^r'); 
legend('Group 1', 'Group 2'); title('C) Principal Coordinates');
xlabel('1st Principal Coordinate'); ylabel('2nd Principal Coordinate');
subplot(2,2,4); 
plot(sum(XYZ(find(shapecluster==1),1:3:51),2), ...
    sum(XYZ(find(shapecluster==1),2:3:51),2), 'ob'); 
grid on; hold on; 
plot(sum(XYZ(find(shapecluster==2),1:3:51),2), ...
    sum(XYZ(find(shapecluster==2),2:3:51),2), '^r');
xlabel('Sum of x-coordinates'); ylabel('Sum of y-coordinates');
legend('Group 1', 'Group 2'); title('D) Sum of Coordinates');

%}

%%
var1 = twist; 
var2 = torglob;
figure; 
plot(var1(cluster_shape == 1), var2(cluster_shape == 1), 'ob'); 
grid on; hold on; 
plot(var1(cluster_shape == 2), var2(cluster_shape == 2), '^r'); 

cluster_var12 = kmeans([var1 var2], 2);
acc = cluster_shape == cluster_var12; acc = sum(acc)/length(acc); 
acc = max(acc, 1-acc); 

title(num2str(acc));

%%
var1 = writhetwist; 
cluster_var1 = kmeans(var1, 2); 
acc = cluster_shape == cluster_var1; acc = sum(acc)/length(acc); 
acc
%acc = max(acc, 1-acc)