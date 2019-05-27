[num, txt] = xlsread('Writhe-pre-post.xlsx');
XYZ = num(:, 6:end); 
N = 32;

shapecluster = num(:,3); 
writhe_preop = num(:,4); 
writhe_postop = num(:,5);
nonsurg = isnan(writhe_postop);

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

abstorsions = abs(torsions); abstorsions2 = abs(torsions2);

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
var1 = topbottomwrithes; 
var2 = twists;
figure; 
plot(var1(shapecluster == 1), var2(shapecluster == 1), 'ob'); 
grid on; hold on; 
plot(var1(shapecluster == 2), var2(shapecluster == 2), '^r'); 