%%{
%{
load('spines_XYZ.mat');

    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    cm = [x;y;z]';
%}

[num, txt] = xlsread('Writhe-pre-post.xlsx');
XYZ = num(:, 6:end); 
shapecluster = num(:,3); 
writhe_preop = num(:,4); 
writhe_postop = num(:,5);
nonsurg = isnan(writhe_postop);
    % get center points 
    x = XYZ(idx, 1:3:51); 
    y = XYZ(idx, 2:3:51); 
    z = XYZ(idx, 3:3:51); 
    cm = [x;y;z]';

%
q = 4;
tau = zeros(1, length(x)-2*q); wr = tau;
d1 = tau; d2 = d1; d3 = d1;
for t = 1:length(tau)
    [tau(t),d,dd,ddd] = lewinerTorsion(cm, t+q, q);
    d1(t) = norm(d); d2(t) = norm(dd); d3(t) = norm(ddd);
    wr(t) = levittWrithe(cm, (0:(2*q)) + t);
end

vertnames = [arrayfun(@(i) ['T' num2str(i)], 1:12, 'UniformOutput', 0), ...
    arrayfun(@(i) ['L' num2str(i)], 1:5, 'UniformOutput', 0)];
vertnames_part = vertnames((q+1):(end-q));

figure('Position', [50 50 1400 700]); 
subplot(1,4,1); plot3(x, y, z, '-o'); grid on; view([0 0]); 
title('sagittal'); xlabel('x'); ylabel('y'); zlabel('z'); 
text(x,y,z,vertnames);
subplot(1,4,2); plot3(x, y, z, '-o'); grid on; view([90 0]);
title('coronal'); xlabel('x'); ylabel('y'); zlabel('z'); 
text(x,y,z,vertnames);
%subplot(1,4,3); plot([d1;d2;d3], z((q+1):(end-q)), '-o'); grid on;
%title('Derivatives'); ylabel('z'); xlabel('derivatives'); ylim([0 1000]);
%legend('1st', '2nd', '3rd');
subplot(1,4,3); plot(d1, z((q+1):(end-q)), '-o'); grid on;
title('1st Derivative'); ylabel('z'); xlabel('2nd derivative'); ylim([0 1000]);
text(d1, z((q+1):(end-q)), vertnames_part);
subplot(1,4,4); plot(d2, z((q+1):(end-q)), '-o'); grid on;
title('2nd Derivative'); ylabel('z'); xlabel('2nd derivative'); ylim([0 1000]);
text(d2, z((q+1):(end-q)), vertnames_part);
%subplot(1,4,4); plot(tau, z((q+1):(end-q)), '-o'); grid on;
%title('Torsion'); ylabel('z'); xlabel('Torsion'); ylim([0 1000]);
%subplot(1,4,4); plot(-wr, z((q+1):(end-q)), '-o'); grid on;
%title('Writhe'); ylabel('z'); xlabel('Writhe'); ylim([0 1000]);
%text(tau, z((q+1):(end-q)), vertnames_part);

%}

%%
%{
load('spines_XYZ.mat');
%Idx = [12, 4, 8, 17, 11, 18, 9];
Idx = [26, 20, 16, 9, 19, 7, 5];
N = length(Idx);

vertnames = [arrayfun(@(i) ['T' num2str(i)], 1:12, 'UniformOutput', 0), ...
    arrayfun(@(i) ['L' num2str(i)], 1:5, 'UniformOutput', 0)];
lnstyle = {'-o', '-*', '-s', '-^', '-p', '-x', '-h', '-d'};

figure('Position', [50 100 1400 600]);
writhe = zeros(1,N); writheabs = zeros(1,N); 
torsions = zeros(1,N); torsionlocs = zeros(1,N); 

for i = 1:N
    idx = Idx(i);
    
    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    cm = [x;y;z]';
    
    writhe(i) = levittWrithe(cm); writheabs(i) = levittWritheAbs(cm);
    
    % get torsion 
    q = 2; % 2 vertebrae above, 2 vertebrae below, current vertebra -> 5 points to fit each cubic
    vertebrae = (1+q):(length(x)-q);
    tau = arrayfun(@(v) lewinerTorsion(cm, v, q), vertebrae);
    [~,tauvert] = max(abs(tau)); torsions(i) = tau(tauvert); % maximum torsion
    torsionlocs(i) = tauvert+q;
    
    subplot(1,3,1); plot3(x, y, z, lnstyle{i}); hold on; 
    subplot(1,3,2); plot3(x, y, z, lnstyle{i}); hold on; 
    subplot(1,3,3); plot3(x, y, z, lnstyle{i}); hold on;
end

lbl = arrayfun(@(i) ['Patient ', num2str(Idx(i)), ': Writhe = ', ...
    num2str(writhe(i)), ', Abs. Writhe = ', ...
    num2str(writheabs(i)), ', Torsion = ', ...
    num2str(torsions(i)), ' at ', vertnames{torsionlocs(i)}], ...
    1:N, 'UniformOutput', 0);
subplot(1,3,1); grid on; view([0 0]); legend(lbl, 'Location', 'southoutside'); 
title('Sagittal'); xlabel('x'); ylabel('y'); zlabel('z'); 
subplot(1,3,2); grid on; view([90 0]); legend(lbl, 'Location', 'southoutside'); 
title('Coronal'); xlabel('x'); ylabel('y'); zlabel('z'); 
subplot(1,3,3); grid on; view([0 90]); legend(lbl, 'Location', 'southoutside'); 
title('Axial'); xlabel('x'); ylabel('y'); zlabel('z'); 
%}