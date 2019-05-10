% Loads N spines from spines_XYZ.mat which contains table spinesXYZ. Gets
% writhe and maximum torsion for each spine. Plots the results for
% clustering. 

load('spines_XYZ.mat');

N = 33;
writhes = zeros(N,1); abswrithes = zeros(N,1); torsions = zeros(N,1); 
for idx = 1:N
    % get center points 
    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    cm = [x;y;z]';
    
    % get writhe 
    writhes(idx) = levittWrithe(cm);
    abswrithes(idx) = levittWritheAbs(cm);
    
    % get torsion 
    q = 2; % 2 vertebrae above, 2 vertebrae below, current vertebra -> 5 points to fit each cubic
    vertebrae = (1+q):(length(x)-q);
    tau = arrayfun(@(v) lewinerTorsion(cm, v, q), vertebrae);
    [~,tauvert] = max(abs(tau)); torsions(idx) = tau(tauvert); % maximum torsion
end

% report magnitude only, not direction 
writhes = abs(writhes); torsions = abs(torsions);

%%
figure; scatter(writhes, torsions, '.k'); grid on;
xlabel('Writhe'); ylabel('Max Torsion'); 

lbl = arrayfun(@(idx) num2str(idx), 1:N, 'UniformOutput', 0);
text(writhes, torsions, lbl);

figure; scatter(writhes, abswrithes, '.k'); grid on;
xlabel('Writhe'); ylabel('Abs Writhe'); 
text(writhes, abswrithes, lbl);

%%
load('spines_XYZ.mat');
N = 33;
XYZ = spinesXYZ{1:end, 2:end};
[C,S,~,~,PE] = pca(XYZ);
PCcluster = kmeans(S(:,1:10), 5);
newcluster = kmeans([writhes, abswrithes, torsions], 5);
fullcluster = kmeans(XYZ, 5);
cluster = PCcluster;

XYZcell = arrayfun(@(idx) [XYZ(idx, 1:3:end); XYZ(idx, 2:3:end); XYZ(idx, 3:3:end)]', ...
    1:N, 'UniformOutput', 0);
avgspine = cell(1, 5); 
figure; 
for k = 1:length(avgspine)
    spn = zeros(size(XYZcell{1}));
    spines = XYZcell(cluster == k);
    for i = 1:length(spines)
        spn = spn + spines{i};
    end
    spn = spn/i;
    for j = 1:3
        subplot(1,3,j); 
        plot3(spn(:,1), spn(:,2), spn(:,3), '-o'); 
        grid on; hold on;
    end
    avgspine{k} = spn;
end
subplot(1,3,1); view([0 0]); 
subplot(1,3,2); view([90 0]); 
subplot(1,3,3); view([0 90]);