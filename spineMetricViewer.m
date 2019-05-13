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

%%
load('spines_XYZ.mat');
Idx = [1,29];
N = length(Idx);

figure('Position', [50 50 1400 700]); 
writhe = zeros(1,N); writheabs = zeros(1,N);

for i = 1:N
    idx = Idx(i);
    
    xyz = spinesXYZ{idx, 2:end};
    x = xyz(1:3:end); 
    y = xyz(2:3:end); 
    z = xyz(3:3:end); 
    cm = [x;y;z]';
    
    writhe(i) = levittWrithe(cm); writheabs(i) = levittWritheAbs(cm);
    subplot(1,2,1); plot3(x, y, z, '-o'); hold on; 
    subplot(1,2,2); plot3(x, y, z, '-o'); hold on; 
end

lbl = arrayfun(@(i) [num2str(Idx(i)), ' | ', num2str(writhe(i)), ' | ', num2str(writheabs(i))], ...
    1:N, 'UniformOutput', 0);
subplot(1,2,1); grid on; view([0 0]); legend(lbl); 
subplot(1,2,2); grid on; view([90 0]); legend(lbl);