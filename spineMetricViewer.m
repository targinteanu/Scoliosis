%{

load('spines_XYZ.mat');

    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    cm = [x;y;z]';

%%
tau = zeros(1, length(x)-2*q); wr = tau;
for t = 1:length(tau)
    tau(t) = lewinerTorsion(cm, t+q, q);
    wr(t) = levittWrithe(cm, (0:(2*q)) + t);
end

figure('Position', [50 50 1400 700]); 
subplot(1,4,1); plot3(x, y, z, '-o'); grid on; view([0 0]);
subplot(1,4,2); plot3(x, y, z, '-o'); grid on; view([90 0]);
subplot(1,4,3); plot(tau, z((q+1):(end-q)), '-o'); grid on;
subplot(1,4,4); plot(-wr, z((q+1):(end-q)), '-o'); grid on;

%}

%%
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